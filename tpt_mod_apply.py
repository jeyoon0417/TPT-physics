#!/usr/bin/env python3
# TPT-physics mod installer: run from the repo root before building.
# Copies the rigid-body module into src/simulation and hooks it into the simulation loop.
import os, re, sys, glob

SRC = os.path.join('src', 'simulation')

FILES = {
'RigidPhysics.h': r'''#pragma once
class Simulation;
// TPT-physics mod: advances rigid bodies made of solid elements. Called once per simulated frame.
void RigidPhysics_Step(Simulation *sim);
''',
'RigidPhysics.cpp': r'''// RigidPhysics.cpp — TPT-physics mod: rigid bodies for solid elements
// Core (namespace rigid) has no TPT dependency; the adapter at the bottom talks to TPT.
#include "RigidPhysics.h"
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <memory>

namespace rigid
{
enum Kind { K_GENERIC = 0, K_GLASS, K_WOOD, K_STONE, K_METAL };
enum Convert { CV_NONE = 0, CV_SHATTER, CV_DUST };

struct Material
{
	bool rigid = false;
	float dens = 1.0f, fr = 0.5f, str = 1.0f;
	int kind = K_GENERIC;
};

struct Env
{
	virtual ~Env() = default;
	virtual bool Blocked(int x, int y) = 0; // static obstacle for rigid particles
	virtual bool WallAt(int x, int y) = 0;  // anchors newly formed bodies
	virtual void Force(float x, float y, float &ax, float &ay) { ax = 0; ay = 0; (void)x; (void)y; }
};

constexpr int PENDING = -2; // drawn, waiting to settle into a body (static)
constexpr int FREE = -1;    // single loose rigid grain (dynamic)

struct Body
{
	std::vector<int> m;
	std::vector<float> qx, qy;
	bool stat = false, alive = true;
	int cool = 8;
};

class World
{
public:
	int W, H, cap;
	float G = 0.10f;       // px / frame^2
	int SUB = 3, ITERS = 2, SETTLE = 12;
	float MAXD = 0.85f, R = 0.49f;

	std::vector<uint8_t> use, seen;
	std::vector<int> type, body, birth, gx, gy;
	std::vector<float> x, y, ox, oy, vx, vy, im, dev, fr, str;
	std::vector<uint8_t> kind;
	std::vector<Body> bodies;
	std::vector<int> tracked;             // ids currently tracked
	std::vector<int> moved;               // output: ids to write back
	std::vector<std::pair<int, int>> convert; // output: (id, Convert)
	int frame = 0;

	World(int w, int h, int c) : W(w), H(h), cap(c)
	{
		use.assign(c, 0); seen.assign(c, 0); type.assign(c, 0); body.assign(c, FREE); birth.assign(c, 0);
		gx.assign(c, 0); gy.assign(c, 0); x.assign(c, 0); y.assign(c, 0); ox.assign(c, 0); oy.assign(c, 0);
		vx.assign(c, 0); vy.assign(c, 0); im.assign(c, 0); dev.assign(c, 0); fr.assign(c, 0.5f); str.assign(c, 1);
		kind.assign(c, 0);
		pixelOwner.assign(size_t(w) * h, -1);
		hashCount.assign(size_t(hw()) * hh(), 0); hashStart.assign(size_t(hw()) * hh() + 1, 0);
	}

	bool Tracked(int i) const { return i >= 0 && i < cap && use[i]; }

	// ---- sync with host ----
	void BeginObserve() { frame++; std::fill(seen.begin(), seen.end(), 0); }

	void Observe(int i, int t, float px, float py, const Material &m)
	{
		if (i < 0 || i >= cap) return;
		seen[i] = 1;
		bool fresh = !use[i] || type[i] != t || std::fabs(px - x[i]) > 1.5f || std::fabs(py - y[i]) > 1.5f;
		if (fresh)
		{
			if (use[i]) Remove(i);
			use[i] = 1; type[i] = t; x[i] = ox[i] = px; y[i] = oy[i] = py; vx[i] = vy[i] = 0;
			body[i] = PENDING; birth[i] = frame; dev[i] = 0; im[i] = 0;
			fr[i] = m.fr; str[i] = m.str; kind[i] = uint8_t(m.kind); dens_[i] = m.dens;
			gx[i] = int(std::floor(px + 0.5f)); gy[i] = int(std::floor(py + 0.5f));
		}
		else { x[i] = px; y[i] = py; }
	}

	// Host adds impulse (e.g. from particle vx/vy set by TPT)
	void Impulse(int i, float ix, float iy)
	{
		if (!Tracked(i) || im[i] == 0) return;
		vx[i] += ix; vy[i] += iy;
	}

	void EndObserve()
	{
		std::vector<int> keep;
		keep.reserve(tracked.size());
		for (int i : tracked) if (use[i] && !seen[i]) Remove(i);
		for (size_t i = 0; i < size_t(cap); i++) if (use[i]) keep.push_back(int(i));
		tracked.swap(keep);
		for (int b : dirty) Resplit(b, nullptr);
		dirty.clear();
	}

	void FormBodies(Env &env)
	{
		std::vector<int> pend;
		for (int i : tracked) if (body[i] == PENDING) pend.push_back(i);
		if (pend.empty()) return;
		for (int i : pend) { int px = rx(i), py = ry(i); if (inGrid(px, py)) pixelOwner[size_t(py) * W + px] = i; }
		std::vector<uint8_t> done(pend.size(), 0);
		std::unordered_map<int, int> idx; idx.reserve(pend.size() * 2);
		for (size_t k = 0; k < pend.size(); k++) idx[pend[k]] = int(k);
		for (size_t k = 0; k < pend.size(); k++)
		{
			if (done[k]) continue;
			std::vector<int> comp, st{ pend[k] }; done[k] = 1;
			bool young = false, anchored = false;
			while (!st.empty())
			{
				int a = st.back(); st.pop_back(); comp.push_back(a);
				if (frame - birth[a] < SETTLE) young = true;
				int px = rx(a), py = ry(a);
				for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++)
					if ((dx || dy) && env.WallAt(px + dx, py + dy)) anchored = true;
				const int nb[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
				for (auto &d : nb)
				{
					int qx = px + d[0], qy = py + d[1];
					if (!inGrid(qx, qy)) continue;
					int o = pixelOwner[size_t(qy) * W + qx];
					if (o < 0) continue;
					auto it = idx.find(o);
					if (it == idx.end() || done[it->second]) continue;
					done[it->second] = 1; st.push_back(o);
				}
			}
			if (young) continue;
			for (int i : comp) { gx[i] = rx(i); gy[i] = ry(i); }
			MakeBody(comp, anchored);
		}
		for (int i : pend) { int px = rx(i), py = ry(i); if (inGrid(px, py)) pixelOwner[size_t(py) * W + px] = -1; }
	}

	void Step(Env &env)
	{
		moved.clear(); convert.clear();
		std::vector<int> dyn;
		for (int i : tracked) if (im[i] > 0) dyn.push_back(i);
		if (dyn.empty()) return;
		const float dt = 1.0f / SUB, maxv = MAXD / dt;
		for (int s = 0; s < SUB; s++)
		{
			for (int i : dyn)
			{
				ox[i] = x[i]; oy[i] = y[i];
				float ax, ay; env.Force(x[i], y[i], ax, ay);
				vx[i] += ax * im[i] * dt; vy[i] += (G + ay * im[i]) * dt;
				float sp = vx[i] * vx[i] + vy[i] * vy[i];
				if (sp > maxv * maxv) { float f = maxv / std::sqrt(sp); vx[i] *= f; vy[i] *= f; }
				x[i] += vx[i] * dt; y[i] += vy[i] * dt;
			}
			for (int i : tracked) if (im[i] == 0) { ox[i] = x[i]; oy[i] = y[i]; }
			BuildHash();
			for (int it = 0; it < ITERS; it++)
			{
				Contacts();
				for (int i : dyn) StaticCollide(env, i);
				ShapeMatch();
			}
			for (int i : dyn) { vx[i] = (x[i] - ox[i]) / dt * 0.999f; vy[i] = (y[i] - oy[i]) / dt * 0.999f; }
		}
		Fractures();
		for (int i : tracked) if (im[i] > 0 || body[i] == FREE) moved.push_back(i);
	}

private:
	std::vector<float> dens_ = std::vector<float>(size_t(cap), 1.0f);
	std::vector<int> pixelOwner, dirty;
	std::vector<int> hashCount, hashStart, sorted;
	static constexpr int HC = 2;
	int hw() const { return W / HC + 2; }
	int hh() const { return H / HC + 2; }
	bool inGrid(int px, int py) const { return px >= 0 && py >= 0 && px < W && py < H; }
	int rx(int i) const { return int(std::floor(x[i] + 0.5f)); }
	int ry(int i) const { return int(std::floor(y[i] + 0.5f)); }

	void Remove(int i)
	{
		if (body[i] >= 0) dirty.push_back(body[i]);
		use[i] = 0; body[i] = FREE; im[i] = 0;
	}

	int MakeBody(const std::vector<int> &mem, bool anchored)
	{
		if (mem.size() == 1)
		{
			int i = mem[0];
			body[i] = FREE; im[i] = anchored ? 0 : 1.0f / dens_[i];
			if (anchored) body[i] = PENDING, birth[i] = -1000000; // stays static
			return -1;
		}
		float M = 0, cx = 0, cy = 0;
		for (int i : mem) { float m = dens_[i]; M += m; cx += gx[i] * m; cy += gy[i] * m; }
		cx /= M; cy /= M;
		Body B;
		B.stat = anchored;
		B.m = mem;
		B.qx.resize(mem.size()); B.qy.resize(mem.size());
		int id = int(bodies.size());
		for (size_t k = 0; k < mem.size(); k++)
		{
			int i = mem[k];
			B.qx[k] = gx[i] - cx; B.qy[k] = gy[i] - cy;
			body[i] = id; im[i] = anchored ? 0 : 1.0f / dens_[i];
			if (anchored) vx[i] = vy[i] = 0;
		}
		// reuse dead slots to keep the vector small
		for (size_t b = 0; b < bodies.size(); b++)
			if (!bodies[b].alive)
			{
				for (int i : mem) body[i] = int(b);
				bodies[b] = std::move(B);
				return int(b);
			}
		bodies.push_back(std::move(B));
		return id;
	}

	std::vector<std::vector<int>> Components(const std::vector<int> &mem, const std::vector<int> *key)
	{
		std::unordered_map<int64_t, int> map; map.reserve(mem.size() * 2);
		std::unordered_map<int, int> kidx; kidx.reserve(mem.size() * 2);
		for (size_t k = 0; k < mem.size(); k++) { map[int64_t(gx[mem[k]]) * 65536 + gy[mem[k]]] = mem[k]; kidx[mem[k]] = int(k); }
		std::unordered_map<int, uint8_t> seenc; seenc.reserve(mem.size() * 2);
		std::vector<std::vector<int>> out;
		for (int s0 : mem)
		{
			if (seenc.count(s0)) continue;
			std::vector<int> comp, st{ s0 }; seenc[s0] = 1;
			while (!st.empty())
			{
				int a = st.back(); st.pop_back(); comp.push_back(a);
				int64_t k = int64_t(gx[a]) * 65536 + gy[a];
				const int64_t d4[4] = { 65536, -65536, 1, -1 };
				for (auto d : d4)
				{
					auto it = map.find(k + d);
					if (it == map.end()) continue;
					int b = it->second;
					if (seenc.count(b)) continue;
					if (key && (*key)[kidx[a]] != (*key)[kidx[b]]) continue;
					seenc[b] = 1; st.push_back(b);
				}
			}
			out.push_back(std::move(comp));
		}
		return out;
	}

	void Resplit(int b, const std::vector<uint8_t> *detachMask)
	{
		if (b < 0 || b >= int(bodies.size()) || !bodies[b].alive) return;
		Body &B = bodies[b];
		B.alive = false;
		bool stat = B.stat;
		std::vector<int> keep;
		for (size_t k = 0; k < B.m.size(); k++)
		{
			int i = B.m[k];
			if (!use[i] || body[i] != b) continue;
			if (detachMask && (*detachMask)[k]) { body[i] = FREE; im[i] = 1.0f / dens_[i]; continue; }
			keep.push_back(i);
		}
		for (auto &c : Components(keep, nullptr)) MakeBody(c, stat);
	}

	void BuildHash()
	{
		int HW = hw(), HH = hh();
		std::fill(hashCount.begin(), hashCount.end(), 0);
		if (sorted.size() < tracked.size()) sorted.resize(tracked.size());
		auto cellOf = [&](int i) {
			int cx = std::clamp(int(x[i] / HC), 0, HW - 1), cy = std::clamp(int(y[i] / HC), 0, HH - 1);
			return cy * HW + cx;
		};
		for (int i : tracked) hashCount[cellOf(i)]++;
		int s = 0;
		for (int c = 0; c < HW * HH; c++) { hashStart[c] = s; s += hashCount[c]; }
		hashStart[HW * HH] = s;
		std::fill(hashCount.begin(), hashCount.end(), 0);
		for (int i : tracked) { int c = cellOf(i); sorted[hashStart[c] + hashCount[c]++] = i; }
	}

	void Contacts()
	{
		int HW = hw(), HH = hh();
		const float D = 2 * R;
		for (int i : tracked)
		{
			int cx = int(x[i] / HC), cy = int(y[i] / HC);
			for (int yy = cy - 1; yy <= cy + 1; yy++)
			{
				if (yy < 0 || yy >= HH) continue;
				for (int xx = cx - 1; xx <= cx + 1; xx++)
				{
					if (xx < 0 || xx >= HW) continue;
					int c = yy * HW + xx;
					for (int s = hashStart[c], e = hashStart[c] + hashCount[c]; s < e; s++)
					{
						int j = sorted[s];
						if (j <= i) continue;
						if (body[i] >= 0 && body[i] == body[j]) continue;
						float ws = im[i] + im[j];
						if (ws <= 0) continue;
						float dx = x[i] - x[j], dy = y[i] - y[j], d2 = dx * dx + dy * dy;
						if (d2 >= D * D) continue;
						float d = std::sqrt(d2);
						if (d < 1e-6f) { dx = 0; dy = -1; d = 1; }
						float nx = dx / d, ny = dy / d, pen = D - d;
						float wi = im[i] / ws, wj = im[j] / ws;
						x[i] += nx * pen * wi; y[i] += ny * pen * wi;
						x[j] -= nx * pen * wj; y[j] -= ny * pen * wj;
						float mu = std::sqrt(fr[i] * fr[j]);
						float rxv = (x[i] - ox[i]) - (x[j] - ox[j]), ryv = (y[i] - oy[i]) - (y[j] - oy[j]);
						float rn = rxv * nx + ryv * ny;
						float tx = rxv - rn * nx, ty = ryv - rn * ny, tl = std::sqrt(tx * tx + ty * ty);
						if (tl > 1e-7f)
						{
							float k = tl < mu * pen * 1.5f ? 1.0f : std::min(1.0f, mu * pen / tl);
							x[i] -= tx * k * wi; y[i] -= ty * k * wi;
							x[j] += tx * k * wj; y[j] += ty * k * wj;
						}
					}
				}
			}
		}
	}

	void StaticCollide(Env &env, int i)
	{
		int cx = rx(i), cy = ry(i);
		for (int yy = cy - 1; yy <= cy + 1; yy++)
			for (int xx = cx - 1; xx <= cx + 1; xx++)
			{
				if (!env.Blocked(xx, yy)) continue;
				float qx = std::clamp(x[i], xx - 0.5f, xx + 0.5f), qy = std::clamp(y[i], yy - 0.5f, yy + 0.5f);
				float dx = x[i] - qx, dy = y[i] - qy, d2 = dx * dx + dy * dy;
				if (d2 >= R * R) continue;
				float nx, ny, pen;
				if (d2 < 1e-10f)
				{
					// centre inside the square: push out along the shallowest axis
					float l = x[i] - (xx - 0.5f), r = (xx + 0.5f) - x[i], t = y[i] - (yy - 0.5f), b = (yy + 0.5f) - y[i];
					float mn = std::min(std::min(l, r), std::min(t, b));
					if (mn == t) { nx = 0; ny = -1; } else if (mn == b) { nx = 0; ny = 1; }
					else if (mn == l) { nx = -1; ny = 0; } else { nx = 1; ny = 0; }
					pen = mn + R;
				}
				else { float d = std::sqrt(d2); nx = dx / d; ny = dy / d; pen = R - d; }
				x[i] += nx * pen; y[i] += ny * pen;
				float rxv = x[i] - ox[i], ryv = y[i] - oy[i], rn = rxv * nx + ryv * ny;
				float tx = rxv - rn * nx, ty = ryv - rn * ny;
				x[i] -= tx * fr[i]; y[i] -= ty * fr[i];
			}
	}

	void ShapeMatch()
	{
		for (auto &B : bodies)
		{
			if (!B.alive || B.stat) continue;
			size_t L = B.m.size();
			float M = 0, cx = 0, cy = 0;
			for (size_t k = 0; k < L; k++) { int i = B.m[k]; float w = 1.0f / im[i]; M += w; cx += x[i] * w; cy += y[i] * w; }
			cx /= M; cy /= M;
			float a00 = 0, a01 = 0, a10 = 0, a11 = 0;
			for (size_t k = 0; k < L; k++)
			{
				int i = B.m[k]; float w = 1.0f / im[i], px = x[i] - cx, py = y[i] - cy;
				a00 += w * px * B.qx[k]; a01 += w * px * B.qy[k]; a10 += w * py * B.qx[k]; a11 += w * py * B.qy[k];
			}
			float th = std::atan2(a10 - a01, a00 + a11), co = std::cos(th), si = std::sin(th);
			for (size_t k = 0; k < L; k++)
			{
				int i = B.m[k];
				float gxp = cx + co * B.qx[k] - si * B.qy[k], gyp = cy + si * B.qx[k] + co * B.qy[k];
				float ddx = x[i] - gxp, ddy = y[i] - gyp, dd = ddx * ddx + ddy * ddy;
				if (dd > dev[i]) dev[i] = dd;
				x[i] = gxp; y[i] = gyp;
			}
		}
	}

	void Crack(int b, int i0, float ratio)
	{
		Body &B = bodies[b];
		size_t L = B.m.size();
		if (L < 4) return;
		int kd = kind[i0];
		float mx = 0, my = 0;
		for (int i : B.m) { mx += gx[i]; my += gy[i]; }
		mx /= L; my /= L;
		float px = float(gx[i0]), py = float(gy[i0]);
		float toC = std::atan2(my - py, mx - px);
		auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return float(seed >> 8) / 16777216.0f; };
		std::vector<float> dirs;
		if (kd == K_GLASS)
		{
			int nl = std::max(1, std::min(5, int(std::ceil(ratio))));
			for (int t = 0; t < nl; t++) dirs.push_back(toC + t * 3.14159265f / nl + (rnd() - 0.5f) * 0.4f);
		}
		else if (kd == K_WOOD)
		{
			float sxx = 0, syy = 0, sxy = 0;
			for (int i : B.m) { float dx = gx[i] - mx, dy = gy[i] - my; sxx += dx * dx; syy += dy * dy; sxy += dx * dy; }
			float g = 0.5f * std::atan2(2 * sxy, sxx - syy);
			px = (px + mx) / 2 + 0.5f; py = (py + my) / 2 + 0.5f;
			dirs.push_back(g + (rnd() - 0.5f) * 0.2f);
			if (ratio > 6) dirs.push_back(g + 1.5708f + (rnd() - 0.5f) * 0.3f);
		}
		else
		{
			dirs.push_back(toC + (rnd() - 0.5f) * 0.8f);
			if (kd == K_STONE && ratio > 4) dirs.push_back(toC + 1.5708f + (rnd() - 0.5f) * 0.8f);
		}
		std::vector<int> key(L, 0);
		std::vector<uint8_t> dust(L, 0);
		for (size_t k = 0; k < L; k++)
		{
			int i = B.m[k], kk = 0;
			for (size_t t = 0; t < dirs.size(); t++)
			{
				float c = std::cos(dirs[t]), s = std::sin(dirs[t]);
				float cr = (gx[i] - px) * s - (gy[i] - py) * c;
				if (cr > 0) kk |= 1 << t;
				if (kd == K_STONE && std::fabs(cr) < 0.6f && rnd() < 0.6f) dust[k] = 1;
			}
			key[k] = kk;
		}
		bool stat = B.stat;
		B.alive = false;
		std::vector<int> keep, keepKey;
		for (size_t k = 0; k < L; k++)
		{
			int i = B.m[k];
			if (!use[i]) continue;
			if (dust[k]) { body[i] = FREE; im[i] = 1.0f / dens_[i]; convert.push_back({ i, CV_DUST }); continue; }
			keep.push_back(i); keepKey.push_back(key[k]);
		}
		for (auto &c : Components(keep, &keepKey))
		{
			if (c.size() <= 2)
			{
				for (int i : c)
				{
					body[i] = FREE; im[i] = 1.0f / dens_[i];
					if (kd == K_GLASS) convert.push_back({ i, CV_SHATTER });
				}
			}
			else
			{
				int nb = MakeBody(c, stat);
				if (nb >= 0) bodies[nb].cool = kd == K_GLASS ? 2 : 10;
			}
		}
	}

	void Fractures()
	{
		int budget = 4;
		for (size_t b = 0; b < bodies.size() && budget > 0; b++)
		{
			Body &B = bodies[b];
			if (!B.alive || B.stat) continue;
			if (B.cool > 0) { B.cool--; continue; }
			int best = -1; float bv = 0;
			for (int i : B.m)
			{
				float s = str[i], r = dev[i] / (s * s);
				if (r > 1 && r > bv) { bv = r; best = i; }
			}
			if (best >= 0) { budget--; Crack(int(b), best, bv); }
		}
		for (int i : tracked) dev[i] = 0;
	}

	uint32_t seed = 12345u;
};
} // namespace rigid

#ifndef RIGID_STANDALONE
// ======================= TPT adapter =======================
#include "Simulation.h"
#include "ElementClasses.h"
#include "SimulationData.h"

namespace
{
using namespace rigid;

struct TptEnv : public Env
{
	Simulation *sim;
	World *world;
	TptEnv(Simulation *s, World *w) : sim(s), world(w) {}

	bool Wall(int x, int y)
	{
		int b = sim->bmap[y / CELL][x / CELL];
		return b && b != WL_FAN && b != WL_STREAM;
	}
	bool Blocked(int x, int y) override
	{
		if (x < CELL || y < CELL || x >= XRES - CELL || y >= YRES - CELL) return true;
		if (Wall(x, y)) return true;
		int r = sim->pmap[y][x];
		if (!r) return false;
		if (world->Tracked(ID(r))) return false;
		auto &elements = SimulationData::CRef().elements;
		return (elements[TYP(r)].Properties & (TYPE_SOLID | TYPE_PART)) != 0;
	}
	bool WallAt(int x, int y) override
	{
		if (x < 0 || y < 0 || x >= XRES || y >= YRES) return false;
		return Wall(x, y);
	}
	void Force(float x, float y, float &ax, float &ay) override
	{
		// air pressure gradient pushes bodies (explosions)
		int cx = int(x + 0.5f) / CELL, cy = int(y + 0.5f) / CELL;
		ax = ay = 0;
		if (cx < 1 || cy < 1 || cx >= XCELLS - 1 || cy >= YCELLS - 1) return;
		const float k = 0.02f;
		ax = (sim->pv[cy][cx - 1] - sim->pv[cy][cx + 1]) * k;
		ay = (sim->pv[cy - 1][cx] - sim->pv[cy + 1][cx]) * k;
	}
};

struct State
{
	std::unique_ptr<World> world;
	Material mat[PT_NUM];
};

State &GetState(Simulation *sim)
{
	static std::unordered_map<Simulation *, std::unique_ptr<State>> states;
	auto &s = states[sim];
	if (!s)
	{
		s = std::make_unique<State>();
		s->world = std::make_unique<World>(XRES, YRES, NPART);
		auto set = [&](int t, float dens, float fr, float str, int kind) {
			s->mat[t].rigid = true; s->mat[t].dens = dens; s->mat[t].fr = fr; s->mat[t].str = str; s->mat[t].kind = kind;
		};
		set(PT_IRON, 7.8f, 0.45f, 1.30f, K_METAL);
		set(PT_METL, 7.8f, 0.45f, 1.20f, K_METAL);
		set(PT_GOLD, 19.3f, 0.40f, 1.00f, K_METAL);
		set(PT_WOOD, 0.6f, 0.65f, 0.95f, K_WOOD);
		set(PT_GLAS, 2.5f, 0.25f, 0.32f, K_GLASS);
		set(PT_BRCK, 2.0f, 0.75f, 0.55f, K_STONE);
	}
	return *s;
}

void Relocate(Simulation *sim, int j, int x0, int y0)
{
	for (int r = 1; r <= 4; r++)
		for (int dy = -r; dy <= r; dy++)
			for (int dx = -r; dx <= r; dx++)
			{
				if (std::max(std::abs(dx), std::abs(dy)) != r) continue;
				int x = x0 + dx, y = y0 + dy;
				if (x < CELL || y < CELL || x >= XRES - CELL || y >= YRES - CELL) continue;
				if (sim->pmap[y][x] || sim->bmap[y / CELL][x / CELL]) continue;
				sim->parts[j].x = float(x); sim->parts[j].y = float(y);
				sim->pmap[y][x] = PMAP(j, sim->parts[j].type);
				return;
			}
}
} // namespace

void RigidPhysics_Step(Simulation *sim)
{
	State &st = GetState(sim);
	World &w = *st.world;
	auto &elements = SimulationData::CRef().elements;

	w.BeginObserve();
	int active = sim->parts.active;
	for (int i = 0; i < active; i++)
	{
		int t = sim->parts[i].type;
		if (t <= 0 || t >= PT_NUM || !st.mat[t].rigid) continue;
		w.Observe(i, t, sim->parts[i].x, sim->parts[i].y, st.mat[t]);
		// velocity given by TPT (e.g. from other elements) becomes an impulse, TPT itself must not move rigid parts
		if (sim->parts[i].vx != 0 || sim->parts[i].vy != 0)
		{
			w.Impulse(i, sim->parts[i].vx, sim->parts[i].vy);
			sim->parts[i].vx = 0; sim->parts[i].vy = 0;
		}
	}
	w.EndObserve();

	TptEnv env(sim, &w);
	w.FormBodies(env);
	w.Step(env);

	// write back positions + pmap
	for (int i : w.moved)
	{
		auto &p = sim->parts[i];
		int oxp = int(p.x + 0.5f), oyp = int(p.y + 0.5f);
		if (oxp >= 0 && oyp >= 0 && oxp < XRES && oyp < YRES && ID(sim->pmap[oyp][oxp]) == i && sim->pmap[oyp][oxp])
			sim->pmap[oyp][oxp] = 0;
	}
	for (int i : w.moved)
	{
		auto &p = sim->parts[i];
		p.x = w.x[i]; p.y = w.y[i]; p.vx = 0; p.vy = 0;
		int nx = int(p.x + 0.5f), ny = int(p.y + 0.5f);
		if (nx < 0 || ny < 0 || nx >= XRES || ny >= YRES) continue;
		int r = sim->pmap[ny][nx];
		if (r && ID(r) != i && !w.Tracked(ID(r)))
		{
			int props = elements[TYP(r)].Properties;
			if (props & (TYPE_LIQUID | TYPE_GAS)) Relocate(sim, ID(r), nx, ny);
		}
		sim->pmap[ny][nx] = PMAP(i, p.type);
	}
	for (auto &c : w.convert)
	{
		int i = c.first;
		auto &p = sim->parts[i];
		int px = int(p.x + 0.5f), py = int(p.y + 0.5f);
		if (c.second == CV_SHATTER && p.type == PT_GLAS) sim->part_change_type(i, px, py, PT_BGLA);
		else if (c.second == CV_DUST && p.type == PT_BRCK) sim->part_change_type(i, px, py, PT_STNE);
	}
}
#endif
''',
}

def fail(msg, path=None):
    print('::error::' + msg)
    if path and os.path.exists(path):
        print('----- ' + path + ' (head) -----')
        print(open(path, encoding='utf-8').read()[:4000])
    sys.exit(1)

# 1) module files
for name, text in FILES.items():
    with open(os.path.join(SRC, name), 'w', encoding='utf-8', newline='\n') as f:
        f.write(text)
print('wrote module files')

# 2) meson: add RigidPhysics.cpp next to Simulation.cpp
mpath = os.path.join(SRC, 'meson.build')
m = open(mpath, encoding='utf-8').read()
if 'RigidPhysics.cpp' not in m:
    mm = re.search(r"^([ \t]*)'Simulation\.cpp',[ \t]*$", m, re.M)
    if not mm:
        fail("could not find 'Simulation.cpp', in " + mpath, mpath)
    indent = mm.group(1)
    m = m[:mm.start()] + indent + "'RigidPhysics.cpp',\n" + m[mm.start():]
    open(mpath, 'w', encoding='utf-8', newline='\n').write(m)
print('patched ' + mpath)

# 3) Simulation.cpp: include + call at the start of Simulation::AfterSim()
spath = os.path.join(SRC, 'Simulation.cpp')
s = open(spath, encoding='utf-8').read()
if 'RigidPhysics_Step' not in s:
    inc = re.search(r'^#include "Simulation\.h"[ \t]*$', s, re.M)
    if not inc:
        fail('could not find #include "Simulation.h" in ' + spath)
    s = s[:inc.end()] + '\n#include "RigidPhysics.h"' + s[inc.end():]
    fn = re.search(r'void\s+Simulation::AfterSim\s*\(\s*\)\s*\{', s)
    if not fn:
        print('Simulation:: functions found:')
        for line in s.splitlines():
            if re.match(r'\s*\w[\w\s\*&:<>]*Simulation::\w+\s*\(', line):
                print('   ' + line.strip())
        fail('could not find Simulation::AfterSim() in ' + spath)
    s = s[:fn.end()] + '\n\tRigidPhysics_Step(this);' + s[fn.end():]
    open(spath, 'w', encoding='utf-8', newline='\n').write(s)
print('patched ' + spath)

# diagnostics for the next steps of the mod
print('----- AfterSim callers -----')
for p in glob.glob('src/**/*.cpp', recursive=True):
    for n, line in enumerate(open(p, encoding='utf-8', errors='replace'), 1):
        if 'AfterSim(' in line or 'BeforeSim(' in line:
            print('%s:%d: %s' % (p, n, line.rstrip()))
print('TPT-physics mod applied OK')

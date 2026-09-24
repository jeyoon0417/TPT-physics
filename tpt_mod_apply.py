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
// Bodies are solved as real rigid bodies (position + angle, mass + inertia, XPBD-style
// position corrections). Only the outline pixels of a body take part in collisions.
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
	float dens = 1.0f, fr = 0.5f, str = 10.0f; // str: survivable impact speed (px/frame)
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
	std::vector<float> qx, qy;          // local offsets = lattice coords - lattice centre (gcx, gcy)
	float gcx = 0, gcy = 0;
	std::unordered_map<int64_t, int> occ; // lattice cell -> member index
	float cx = 0, cy = 0, th = 0, co = 1, si = 0;
	float vx = 0, vy = 0, w = 0;
	float pcx = 0, pcy = 0, pth = 0;
	float pvx = 0, pvy = 0, pw = 0;
	float M = 1, I = 1, rmax = 1;
	float fax = 0, fay = 0, torque = 0; // external force per frame (per unit mass / inertia)
	bool stat = false, alive = true, asleep = false;
	int cool = 8, still = 0;
	float impact = 0, recent = 0;
	float adx = 0, ady = 0, adth = 0; int an = 0; // correction accumulator (Jacobi)
	float avx = 0, avy = 0, aw = 0; int avn = 0;  // velocity impulse accumulator
	void SetAngle(float a) { th = a; co = std::cos(a); si = std::sin(a); }
};

class World
{
public:
	int W, H, cap;
	float G = 0.20f;         // px / frame^2
	int SUB = 3, ITERS = 16, SETTLE = 12, VPASS = 10;
	float MAXD = 0.8f, R = 0.49f, VMAX = 12.0f;

	std::vector<uint8_t> use, seen, boundary;
	std::vector<int> type, body, bidx, birth, gx, gy;
	std::vector<float> x, y, ox, oy, vx, vy, dens, fr, str, corr;
	std::vector<uint8_t> kind;
	std::vector<Body> bodies;
	std::vector<int> tracked;
	std::vector<int> moved;
	std::vector<std::pair<int, int>> convert;
	int frame = 0;

	World(int w, int h, int c) : W(w), H(h), cap(c)
	{
		use.assign(c, 0); seen.assign(c, 0); boundary.assign(c, 0);
		type.assign(c, 0); body.assign(c, FREE); bidx.assign(c, 0); birth.assign(c, 0); gx.assign(c, 0); gy.assign(c, 0);
		x.assign(c, 0); y.assign(c, 0); ox.assign(c, 0); oy.assign(c, 0); vx.assign(c, 0); vy.assign(c, 0);
		dens.assign(c, 1); fr.assign(c, 0.5f); str.assign(c, 10); corr.assign(c, 0); kind.assign(c, 0);
		pixelOwner.assign(size_t(w) * h, -1);
		fdx.assign(c, 0); fdy.assign(c, 0); fn.assign(c, 0);
		pendPix.assign(size_t(w) * h, 0);
		hashCount.assign(size_t(hw()) * hh(), 0); hashStart.assign(size_t(hw()) * hh() + 1, 0);
	}

	bool Tracked(int i) const { return i >= 0 && i < cap && use[i]; }

	// ---------------- sync with host ----------------
	void BeginObserve() { frame++; std::fill(seen.begin(), seen.end(), 0); }

	void Observe(int i, int t, float px, float py, const Material &m)
	{
		if (i < 0 || i >= cap) return;
		seen[i] = 1;
		bool fresh = !use[i] || type[i] != t || std::fabs(px - x[i]) > 1.5f || std::fabs(py - y[i]) > 1.5f;
		if (!fresh) return;
		if (use[i]) Remove(i);
		use[i] = 1; type[i] = t; x[i] = ox[i] = px; y[i] = oy[i] = py; vx[i] = vy[i] = 0;
		body[i] = PENDING; birth[i] = frame; corr[i] = 0;
		dens[i] = m.dens; fr[i] = m.fr; str[i] = m.str; kind[i] = uint8_t(m.kind);
		gx[i] = int(std::floor(px + 0.5f)); gy[i] = int(std::floor(py + 0.5f));
	}

	void Impulse(int i, float ix, float iy)
	{
		if (!Tracked(i)) return;
		int b = body[i];
		if (b >= 0)
		{
			Body &B = bodies[b];
			if (B.stat) return;
			if (B.asleep) Wake(b);
			float m = dens[i];
			B.vx += ix * m / B.M; B.vy += iy * m / B.M;
			B.w += ((x[i] - B.cx) * iy - (y[i] - B.cy) * ix) * m / B.I;
		}
		else if (b == FREE) { vx[i] += ix; vy[i] += iy; }
	}

	void EndObserve()
	{
		for (int i : tracked) if (use[i] && !seen[i]) Remove(i);
		std::vector<int> keep;
		keep.reserve(tracked.size() + 64);
		for (size_t i = 0; i < size_t(cap); i++) if (use[i]) keep.push_back(int(i));
		tracked.swap(keep);
		std::sort(dirty.begin(), dirty.end());
		dirty.erase(std::unique(dirty.begin(), dirty.end()), dirty.end());
		for (int b : dirty) Split(b, nullptr, nullptr, 0);
		dirty.clear();
	}

	void FormBodies(Env &env)
	{
		std::vector<int> pend;
		for (int i : tracked) if (body[i] == PENDING) pend.push_back(i);
		if (pend.empty()) return;
		for (int i : pend) if (inGrid(gx[i], gy[i])) pixelOwner[size_t(gy[i]) * W + gx[i]] = i;
		std::unordered_map<int, int> idx; idx.reserve(pend.size() * 2);
		for (size_t k = 0; k < pend.size(); k++) idx[pend[k]] = int(k);
		std::vector<uint8_t> done(pend.size(), 0);
		for (size_t k = 0; k < pend.size(); k++)
		{
			if (done[k]) continue;
			std::vector<int> comp, st{ pend[k] }; done[k] = 1;
			bool young = false, anchored = false;
			while (!st.empty())
			{
				int a = st.back(); st.pop_back(); comp.push_back(a);
				if (frame - birth[a] < SETTLE) young = true;
				for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++)
					if ((dx || dy) && env.WallAt(gx[a] + dx, gy[a] + dy)) anchored = true;
				const int nb[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
				for (auto &d : nb)
				{
					int qx = gx[a] + d[0], qy = gy[a] + d[1];
					if (!inGrid(qx, qy)) continue;
					int o = pixelOwner[size_t(qy) * W + qx];
					if (o < 0) continue;
					auto it = idx.find(o);
					if (it == idx.end() || done[it->second]) continue;
					done[it->second] = 1; st.push_back(o);
				}
			}
			if (young) continue;
			std::vector<float> pvx(comp.size(), 0), pvy(comp.size(), 0);
			MakeBody(comp, anchored, pvx, pvy);
		}
		for (int i : pend) if (inGrid(gx[i], gy[i])) pixelOwner[size_t(gy[i]) * W + gx[i]] = -1;
	}

	void Step(Env &env)
	{
		moved.clear(); convert.clear();
		for (int b : wakeNext) if (b >= 0 && b < int(bodies.size()) && bodies[b].alive) Wake(b);
		wakeNext.clear();
		if (frame % 30 == 0)
			for (size_t b = 0; b < bodies.size(); b++)
				if (bodies[b].alive && bodies[b].asleep) { Wake(int(b)); bodies[b].still = 13; }

		for (int i : pendList) if (inGrid(gx[i], gy[i])) pendPix[size_t(gy[i]) * W + gx[i]] = 0;
		pendList.clear();
		for (int i : tracked) if (body[i] == PENDING) { pendList.push_back(i); if (inGrid(gx[i], gy[i])) pendPix[size_t(gy[i]) * W + gx[i]] = 1; }
		std::vector<int> dynB, freeP;
		for (size_t b = 0; b < bodies.size(); b++) if (bodies[b].alive && !bodies[b].stat && !bodies[b].asleep) dynB.push_back(int(b));
		for (int i : tracked) if (body[i] == FREE) freeP.push_back(i);
		if (dynB.empty() && freeP.empty()) return;

		// external forces (air pressure), once per frame
		for (int b : dynB)
		{
			Body &B = bodies[b];
			B.fax = B.fay = B.torque = 0;
			for (size_t k = 0; k < B.m.size(); k += 1)
			{
				int i = B.m[k];
				if (!boundary[i]) continue;
				float px, py; Pos(i, px, py);
				float ax, ay; env.Force(px, py, ax, ay);
				B.fax += ax; B.fay += ay;
				B.torque += (px - B.cx) * ay - (py - B.cy) * ax;
			}
			B.fax /= B.M; B.fay /= B.M; B.torque /= B.I;
		}

		float vmax = 0;
		for (int b : dynB) { Body &B = bodies[b]; vmax = std::max(vmax, std::hypot(B.vx, B.vy) + std::fabs(B.w) * B.rmax); }
		for (int i : freeP) vmax = std::max(vmax, std::hypot(vx[i], vy[i]));
		vmax = std::min(vmax + G * 2, VMAX);
		SUB = std::clamp(int(std::ceil(vmax / MAXD)), 2, 16);
		const float dt = 1.0f / SUB, maxv = MAXD / dt;

		for (int s = 0; s < SUB; s++)
		{
			// integrate
			for (int b : dynB)
			{
				Body &B = bodies[b];
				B.pcx = B.cx; B.pcy = B.cy; B.pth = B.th;
				B.vx += B.fax * dt; B.vy += (G + B.fay) * dt; B.w += B.torque * dt;
				float sp = std::hypot(B.vx, B.vy) + std::fabs(B.w) * B.rmax;
				if (sp > maxv) { float f = maxv / sp; B.vx *= f; B.vy *= f; B.w *= f; }
				B.pvx = B.vx; B.pvy = B.vy; B.pw = B.w;
				B.cx += B.vx * dt; B.cy += B.vy * dt; B.SetAngle(B.th + B.w * dt);
			}
			for (int i : freeP)
			{
				vy[i] += G * dt;
				float sp = std::hypot(vx[i], vy[i]);
				if (sp > maxv) { vx[i] *= maxv / sp; vy[i] *= maxv / sp; }
				x[i] += vx[i] * dt; y[i] += vy[i] * dt;
			}
			// collision candidates: outline pixels of every body + loose grains + pending pixels
			hashList.clear();
			for (int i : tracked)
			{
				if (body[i] == PENDING || (body[i] >= 0 && !boundary[i])) continue;
				float px, py; Pos(i, px, py);
				x[i] = px; y[i] = py;
				hashList.push_back(i);
			}
			// start-of-substep positions (for friction)
			for (int b : dynB)
			{
				Body &B = bodies[b];
				float c0 = std::cos(B.pth), s0 = std::sin(B.pth);
				for (size_t k = 0; k < B.m.size(); k++)
				{
					int i = B.m[k];
					if (!boundary[i]) continue;
					ox[i] = B.pcx + c0 * B.qx[k] - s0 * B.qy[k]; oy[i] = B.pcy + s0 * B.qx[k] + c0 * B.qy[k];
				}
			}
			for (int i : freeP) { ox[i] = x[i] - vx[i] * dt; oy[i] = y[i] - vy[i] * dt; }
			BuildHash();
			conts.clear();
			for (int it = 0; it < ITERS; it++)
			{
				recording = it == ITERS - 1;
				Contacts();
				for (int b : dynB)
				{
					Body &B = bodies[b];
					for (int i : B.m) if (boundary[i]) StaticCollide(env, i);
				}
				for (int i : freeP) StaticCollide(env, i);
				Commit(dynB, freeP);
			}
			// velocities from positions
			for (int b : dynB)
			{
				Body &B = bodies[b];
				B.vx = (B.cx - B.pcx) / dt; B.vy = (B.cy - B.pcy) / dt; B.w = (B.th - B.pth) / dt;
			}
			for (int i : freeP) { float px = x[i], py = y[i]; vx[i] = (px - ox[i]) / dt; vy[i] = (py - oy[i]) / dt; }
			recording = false;
			for (int b : dynB)
			{
				Body &B = bodies[b];
				B.impact += std::hypot(B.vx - B.pvx, B.vy - B.pvy) + std::fabs(B.w - B.pw) * B.rmax * 0.5f;
			}
		}
		// final particle positions
		for (int b : dynB)
		{
			Body &B = bodies[b];
			for (int i : B.m) { Pos(i, x[i], y[i]); moved.push_back(i); }
		}
		for (int i : freeP) moved.push_back(i);
		Fractures();
		Sleep();
	}

private:
	std::vector<int> pixelOwner, dirty, wakeNext, hashList;
	std::vector<float> fdx, fdy; std::vector<int> fn;
	std::vector<uint8_t> pendPix; std::vector<int> pendList;
	struct Cont { int i, b, j; float nx, ny, qx, qy; }; // b: body hit (-1 static/none), j: grain hit, q: world point on b
	std::vector<Cont> conts; bool recording = false;
	std::vector<int> hashCount, hashStart, sorted;
	static constexpr int HC = 2;
	uint32_t seed = 12345u;
	int hw() const { return W / HC + 2; }
	int hh() const { return H / HC + 2; }
	bool inGrid(int px, int py) const { return px >= 0 && py >= 0 && px < W && py < H; }
	float rnd() { seed = seed * 1664525u + 1013904223u; return float(seed >> 8) / 16777216.0f; }

	// ---------------- body helpers ----------------
	void Pos(int i, float &px, float &py) const
	{
		int b = body[i];
		if (b >= 0)
		{
			const Body &B = bodies[b];
			int k = bidx[i];
			px = B.cx + B.co * B.qx[k] - B.si * B.qy[k];
			py = B.cy + B.si * B.qx[k] + B.co * B.qy[k];
		}
		else { px = x[i]; py = y[i]; }
	}
	bool Movable(int i) const
	{
		int b = body[i];
		if (b >= 0) { const Body &B = bodies[b]; return !B.stat && !B.asleep; }
		return b == FREE;
	}
	// generalised inverse mass of particle i for a push along n
	float InvMass(int i, float px, float py, float nx, float ny) const
	{
		int b = body[i];
		if (b >= 0)
		{
			const Body &B = bodies[b];
			if (B.stat || B.asleep) return 0;
			float rn = (px - B.cx) * ny - (py - B.cy) * nx;
			return 1.0f / B.M + rn * rn / B.I;
		}
		return b == FREE ? 1.0f / dens[i] : 0.0f;
	}
	// accumulate a positional impulse (dx, dy) = n * lambda at particle i; applied averaged per iteration,
	// so the order in which contacts are visited does not bias the result (no fake torques)
	void Apply(int i, float px, float py, float dx, float dy)
	{
		int b = body[i];
		if (b >= 0)
		{
			Body &B = bodies[b];
			if (B.stat || B.asleep) return;
			float rx = px - B.cx, ry = py - B.cy;
			B.adx += dx / B.M; B.ady += dy / B.M; B.adth += (rx * dy - ry * dx) / B.I; B.an++;
		}
		else if (b == FREE) { fdx[i] += dx; fdy[i] += dy; fn[i]++; }
	}
	void Commit(const std::vector<int> &dynB, const std::vector<int> &freeP)
	{
		for (int b : dynB)
		{
			Body &B = bodies[b];
			if (!B.an) continue;
			float k = 1.0f / B.an;
			B.cx += B.adx * k; B.cy += B.ady * k; B.SetAngle(B.th + B.adth * k);
			B.adx = B.ady = B.adth = 0; B.an = 0;
		}
		for (int i : freeP)
		{
			if (!fn[i]) continue;
			x[i] += fdx[i] / fn[i]; y[i] += fdy[i] / fn[i];
			fdx[i] = fdy[i] = 0; fn[i] = 0;
		}
	}

	void Wake(int b)
	{
		Body &B = bodies[b];
		if (!B.asleep) return;
		B.asleep = false; B.still = 0;
	}
	void Sleep()
	{
		for (auto &B : bodies)
		{
			if (!B.alive || B.stat || B.asleep) continue;
			float sp = std::hypot(B.vx, B.vy) + std::fabs(B.w) * B.rmax;
			if (sp < 0.06f) B.still++; else B.still = 0;
			if (B.still > 15) { B.asleep = true; B.vx = B.vy = B.w = 0; }
		}
	}

	void Remove(int i)
	{
		if (body[i] >= 0) dirty.push_back(body[i]);
		use[i] = 0; body[i] = FREE;
	}

	static int64_t LKey(int ax, int ay) { return int64_t(ax) * 1048576 + ay; }

	// build a body from particles; pvx/pvy = current velocity of each particle (for momentum);
	// theta = orientation of the lattice (0 for freshly drawn pixels, parent's angle for fragments)
	int MakeBody(const std::vector<int> &mem, bool anchored, const std::vector<float> &pvx, const std::vector<float> &pvy, float theta = 0)
	{
		if (mem.size() == 1)
		{
			int i = mem[0];
			if (anchored) { body[i] = PENDING; birth[i] = -1000000; return -1; }
			body[i] = FREE; vx[i] = pvx[0]; vy[i] = pvy[0]; boundary[i] = 1;
			return -1;
		}
		Body B;
		B.stat = anchored;
		B.m = mem;
		float M = 0, cx = 0, cy = 0, gcx = 0, gcy = 0, mvx = 0, mvy = 0;
		for (size_t k = 0; k < mem.size(); k++)
		{
			int i = mem[k]; float m = dens[i];
			M += m; cx += x[i] * m; cy += y[i] * m; gcx += gx[i] * m; gcy += gy[i] * m;
			mvx += pvx[k] * m; mvy += pvy[k] * m;
		}
		cx /= M; cy /= M; gcx /= M; gcy /= M;
		B.M = M; B.cx = cx; B.cy = cy; B.gcx = gcx; B.gcy = gcy; B.SetAngle(theta);
		B.pcx = cx; B.pcy = cy; B.pth = theta;
		B.qx.resize(mem.size()); B.qy.resize(mem.size());
		B.occ.reserve(mem.size() * 2);
		float I = 0, L = 0, rmax = 0;
		for (size_t k = 0; k < mem.size(); k++)
		{
			int i = mem[k]; float m = dens[i];
			B.qx[k] = gx[i] - gcx; B.qy[k] = gy[i] - gcy;
			B.occ[LKey(gx[i], gy[i])] = int(k);
			float rx = B.qx[k], ry = B.qy[k];
			I += m * (rx * rx + ry * ry + 1.0f / 6.0f);
			float wx = x[i] - cx, wy = y[i] - cy;
			L += m * (wx * (pvy[k] - mvy / M) - wy * (pvx[k] - mvx / M));
			rmax = std::max(rmax, std::sqrt(rx * rx + ry * ry));
		}
		B.I = std::max(I, 1e-3f); B.rmax = rmax + 0.5f;
		if (!anchored) { B.vx = mvx / M; B.vy = mvy / M; B.w = L / B.I; }
		int id = -1;
		for (size_t b = 0; b < bodies.size(); b++) if (!bodies[b].alive) { id = int(b); break; }
		if (id < 0) { id = int(bodies.size()); bodies.emplace_back(); }
		for (size_t k = 0; k < mem.size(); k++)
		{
			int i = mem[k];
			body[i] = id; bidx[i] = int(k);
			boundary[i] = !(B.occ.count(LKey(gx[i] + 1, gy[i])) && B.occ.count(LKey(gx[i] - 1, gy[i])) &&
				B.occ.count(LKey(gx[i], gy[i] + 1)) && B.occ.count(LKey(gx[i], gy[i] - 1)));
		}
		bodies[id] = std::move(B);
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

	// split body b into pieces: key (per member) separates pieces, detach (per member) frees particles
	void Split(int b, const std::vector<int> *key, const std::vector<uint8_t> *detach, int kd)
	{
		if (b < 0 || b >= int(bodies.size()) || !bodies[b].alive) return;
		Body &B = bodies[b];
		bool stat = B.stat;
		std::vector<int> keep, keepKey;
		std::unordered_map<int, std::pair<float, float>> vel; vel.reserve(B.m.size() * 2);
		for (size_t k = 0; k < B.m.size(); k++)
		{
			int i = B.m[k];
			if (!use[i] || body[i] != b) continue;
			float px, py; Pos(i, px, py); x[i] = px; y[i] = py;
			float pvx = stat ? 0 : B.vx - B.w * (py - B.cy), pvy = stat ? 0 : B.vy + B.w * (px - B.cx);
			vel[i] = { pvx, pvy };
			if (detach && (*detach)[k])
			{
				body[i] = FREE; boundary[i] = 1; vx[i] = pvx; vy[i] = pvy;
				if (kd == K_STONE) convert.push_back({ i, CV_DUST });
				continue;
			}
			keep.push_back(i); keepKey.push_back(key ? (*key)[k] : 0);
		}
		B.alive = false;
		float theta = B.th;
		for (auto &c : Components(keep, key ? &keepKey : nullptr))
		{
			std::vector<float> cvx(c.size()), cvy(c.size());
			for (size_t k = 0; k < c.size(); k++) { cvx[k] = vel[c[k]].first; cvy[k] = vel[c[k]].second; }
			if (key && c.size() <= 2)
			{
				for (size_t k = 0; k < c.size(); k++)
				{
					int i = c[k]; body[i] = FREE; boundary[i] = 1; vx[i] = cvx[k]; vy[i] = cvy[k];
					if (kd == K_GLASS) convert.push_back({ i, CV_SHATTER });
				}
				continue;
			}
			int nb = MakeBody(c, stat, cvx, cvy, theta);
			if (nb >= 0 && key) bodies[nb].cool = kd == K_GLASS ? 2 : 10;
		}
	}

	// ---------------- collision ----------------
	void BuildHash()
	{
		int HW = hw(), HH = hh();
		std::fill(hashCount.begin(), hashCount.end(), 0);
		if (sorted.size() < hashList.size()) sorted.resize(hashList.size());
		auto cellOf = [&](int i) {
			int cx = std::clamp(int(x[i] / HC), 0, HW - 1), cy = std::clamp(int(y[i] / HC), 0, HH - 1);
			return cy * HW + cx;
		};
		for (int i : hashList) hashCount[cellOf(i)]++;
		int s = 0;
		for (int c = 0; c < HW * HH; c++) { hashStart[c] = s; s += hashCount[c]; }
		hashStart[HW * HH] = s;
		std::fill(hashCount.begin(), hashCount.end(), 0);
		for (int i : hashList) { int c = cellOf(i); sorted[hashStart[c] + hashCount[c]++] = i; }
	}

	// velocity of the material point (px, py) of particle i's owner
	void PointVel(int i, float px, float py, float &ux, float &uy) const
	{
		int b = body[i];
		if (b >= 0) { const Body &B = bodies[b]; ux = B.vx - B.w * (py - B.cy); uy = B.vy + B.w * (px - B.cx); }
		else if (b == FREE) { ux = vx[i]; uy = vy[i]; }
		else { ux = uy = 0; }
	}
	void BodyVel(const Body &B, float px, float py, float &ux, float &uy) const
	{
		if (B.stat || B.asleep) { ux = uy = 0; return; }
		ux = B.vx - B.w * (py - B.cy); uy = B.vy + B.w * (px - B.cx);
	}
	void VelImpulse(int i, float px, float py, float jx, float jy)
	{
		int b = body[i];
		if (b >= 0) { Body &B = bodies[b]; BodyVelImpulse(B, px, py, jx, jy); }
		else if (b == FREE) { vx[i] += jx / dens[i]; vy[i] += jy / dens[i]; }
	}
	void BodyVelImpulse(Body &B, float px, float py, float jx, float jy)
	{
		if (B.stat || B.asleep) return;
		B.vx += jx / B.M; B.vy += jy / B.M; B.w += ((px - B.cx) * jy - (py - B.cy) * jx) / B.I;
	}
	// inelastic contacts (XPBD velocity step): every touching point ends the substep with zero normal
	// velocity, so position corrections never turn into bounces; iterated so the support forces spread
	// correctly over all touching pixels (no fake torque from visiting order). Friction limited by mu * support.
	void VelocityPass(float dt)
	{
		std::vector<float> accN(conts.size(), 0.0f), accT(conts.size(), 0.0f);
		std::vector<int> order(conts.size());
		for (size_t k = 0; k < order.size(); k++) order[k] = int(k);
		for (int pass = 0; pass < VPASS; pass++)
		{
			// alternate / shuffle visiting order so no side of a body is systematically favoured
			for (size_t k = order.size(); k > 1; k--) std::swap(order[k - 1], order[size_t(rnd() * k) % k]);
			for (int ci : order)
			{
				auto &c = conts[ci];
				float px, py; Pos(c.i, px, py);
				float ux, uy; PointVel(c.i, px, py, ux, uy);
				float ox_ = 0, oy_ = 0, wo = 0;
				float qx = c.qx, qy = c.qy;
				if (c.b >= 0) BodyVel(bodies[c.b], qx, qy, ox_, oy_);
				else if (c.j >= 0) { ox_ = vx[c.j]; oy_ = vy[c.j]; qx = x[c.j]; qy = y[c.j]; }
				float rvx = ux - ox_, rvy = uy - oy_;
				float vn = rvx * c.nx + rvy * c.ny;
				float wi = InvMass(c.i, px, py, c.nx, c.ny);
				if (c.b >= 0) wo = BodyInvMass(bodies[c.b], qx, qy, c.nx, c.ny);
				else if (c.j >= 0) wo = 1.0f / dens[c.j];
				if (wi + wo <= 0) continue;
				float P = -vn / (wi + wo);
				accN[ci] += P;
				VelImpulse(c.i, px, py, c.nx * P, c.ny * P);
				if (c.b >= 0) BodyVelImpulse(bodies[c.b], qx, qy, -c.nx * P, -c.ny * P);
				else if (c.j >= 0) { vx[c.j] -= c.nx * P / dens[c.j]; vy[c.j] -= c.ny * P / dens[c.j]; }
				// friction
				rvx += c.nx * P * wi; rvy += c.ny * P * wi; // approx. updated relative velocity
				float vt_x = rvx - (rvx * c.nx + rvy * c.ny) * c.nx, vt_y = rvy - (rvx * c.nx + rvy * c.ny) * c.ny;
				float tl = std::sqrt(vt_x * vt_x + vt_y * vt_y);
				if (tl < 1e-6f) continue;
				float tx = vt_x / tl, ty = vt_y / tl;
				float wti = InvMass(c.i, px, py, tx, ty), wto = 0;
				if (c.b >= 0) wto = BodyInvMass(bodies[c.b], qx, qy, tx, ty);
				else if (c.j >= 0) wto = 1.0f / dens[c.j];
				if (wti + wto <= 0) continue;
				float mu = fr[c.i];
				float support = std::fabs(accN[ci]) + G * dt / (wi + wo) * std::max(0.0f, -c.ny);
				float want = tl / (wti + wto);
				float lim = mu * support;
				float newT = std::min(accT[ci] + want, lim);
				float Pt = newT - accT[ci];
				if (Pt <= 0) continue;
				accT[ci] = newT;
				VelImpulse(c.i, px, py, -tx * Pt, -ty * Pt);
				if (c.b >= 0) BodyVelImpulse(bodies[c.b], qx, qy, tx * Pt, ty * Pt);
				else if (c.j >= 0) { vx[c.j] += tx * Pt / dens[c.j]; vy[c.j] += ty * Pt / dens[c.j]; }
			}
		}
	}

	// generalised inverse mass of body b at world point (px, py) along n
	float BodyInvMass(const Body &B, float px, float py, float nx, float ny) const
	{
		if (B.stat || B.asleep) return 0;
		float rn = (px - B.cx) * ny - (py - B.cy) * nx;
		return 1.0f / B.M + rn * rn / B.I;
	}
	void BodyApply(Body &B, float px, float py, float dx, float dy)
	{
		if (B.stat || B.asleep) return;
		float rx = px - B.cx, ry = py - B.cy;
		B.adx += dx / B.M; B.ady += dy / B.M; B.adth += (rx * dy - ry * dx) / B.I; B.an++;
	}

	// particle i (circle) against the pixels of body b (unit squares in b's lattice frame)
	void CollideWithBody(int i, int b)
	{
		Body &B = bodies[b];
		float px, py; Pos(i, px, py);
		// into b's lattice coordinates
		float dxw = px - B.cx, dyw = py - B.cy;
		float lx = B.co * dxw + B.si * dyw + B.gcx, ly = -B.si * dxw + B.co * dyw + B.gcy;
		int cx = int(std::floor(lx + 0.5f)), cy = int(std::floor(ly + 0.5f));
		for (int yy = cy - 1; yy <= cy + 1; yy++)
			for (int xx = cx - 1; xx <= cx + 1; xx++)
			{
				if (!B.occ.count(LKey(xx, yy))) continue;
				float qx = std::clamp(lx, xx - 0.5f, xx + 0.5f), qy = std::clamp(ly, yy - 0.5f, yy + 0.5f);
				float dx = lx - qx, dy = ly - qy, d2 = dx * dx + dy * dy;
				float nlx, nly, pen;
				if (d2 < 1e-12f)
				{
					// centre inside b's pixel: leave through the free side facing where the particle came from
					float odx = ox[i] - B.cx, ody = oy[i] - B.cy;
					float l0x = B.co * odx + B.si * ody + B.gcx, l0y = -B.si * odx + B.co * ody + B.gcy;
					float best = 1e9f, bestDist = 0; nlx = 0; nly = -1; pen = 0;
					const int dirs[4][2] = { {0,-1},{0,1},{-1,0},{1,0} };
					for (auto &d : dirs)
					{
						if (B.occ.count(LKey(xx + d[0], yy + d[1]))) continue;
						float dist = d[0] ? (d[0] > 0 ? xx + 0.5f - lx : lx - (xx - 0.5f)) : (d[1] > 0 ? yy + 0.5f - ly : ly - (yy - 0.5f));
						float ex = xx + d[0] * (0.5f + R), ey = yy + d[1] * (0.5f + R);
						if (!d[0]) ex = lx;
						if (!d[1]) ey = ly;
						float score = (ex - l0x) * (ex - l0x) + (ey - l0y) * (ey - l0y);
						if (score < best) { best = score; bestDist = dist; nlx = float(d[0]); nly = float(d[1]); }
					}
					if (best >= 1e9f) continue;
					pen = bestDist + R;
				}
				else
				{
					if (d2 >= R * R) continue;
					// closest point is a corner shared with a neighbouring pixel: the surface is flat there, skip
					bool cornerX = qx != lx, cornerY = qy != ly;
					if (cornerX && cornerY)
					{
						int sx = lx > qx ? 1 : -1, sy = ly > qy ? 1 : -1;
						if (B.occ.count(LKey(xx + sx, yy)) || B.occ.count(LKey(xx, yy + sy))) continue;
					}
					float d = std::sqrt(d2); nlx = dx / d; nly = dy / d; pen = R - d;
				}
				// back to world
				float nx = B.co * nlx - B.si * nly, ny = B.si * nlx + B.co * nly;
				float qwx = B.cx + B.co * (qx - B.gcx) - B.si * (qy - B.gcy);
				float qwy = B.cy + B.si * (qx - B.gcx) + B.co * (qy - B.gcy);
				float wi = InvMass(i, px, py, nx, ny), wb = BodyInvMass(B, qwx, qwy, nx, ny);
				if (wi + wb <= 0) continue;
				if (B.asleep && pen > 0.05f) wakeNext.push_back(b);
				float lam = pen / (wi + wb);
				Apply(i, px, py, nx * lam, ny * lam);
				BodyApply(B, qwx, qwy, -nx * lam, -ny * lam);
				if (recording) conts.push_back({ i, b, -1, nx, ny, qwx, qwy });
				corr[i] = std::max(corr[i], lam * wi);
				int k = B.occ[LKey(xx, yy)]; corr[B.m[k]] = std::max(corr[B.m[k]], lam * wb);
				// friction: relative sliding of the two surfaces during this substep
				float c0 = std::cos(B.pth), s0 = std::sin(B.pth);
				float lqx = qx - B.gcx, lqy = qy - B.gcy;
				float oqx = B.pcx + c0 * lqx - s0 * lqy, oqy = B.pcy + s0 * lqx + c0 * lqy;
				float rxv = (px - ox[i]) - (qwx - oqx), ryv = (py - oy[i]) - (qwy - oqy);
				float rn = rxv * nx + ryv * ny;
				float tx = rxv - rn * nx, ty = ryv - rn * ny, tl = std::sqrt(tx * tx + ty * ty);
				if (tl < 1e-7f) continue;
				tx /= tl; ty /= tl;
				float wti = InvMass(i, px, py, tx, ty), wtb = BodyInvMass(B, qwx, qwy, tx, ty);
				if (wti + wtb <= 0) continue;
				float mu = std::sqrt(fr[i] * fr[B.m[k]]);
				float lt = tl / (wti + wtb);
				if (tl > mu * pen * 1.5f) lt = std::min(lt, mu * lam);
				Apply(i, px, py, -tx * lt, -ty * lt);
				BodyApply(B, qwx, qwy, tx * lt, ty * lt);
			}
	}

	void Contacts()
	{
		int HW = hw(), HH = hh();
		const float D = 2 * R;
		std::vector<int> near;
		for (int i : hashList)
		{
			if (!Movable(i)) continue;
			int cx = std::clamp(int(x[i] / HC), 0, HW - 1), cy = std::clamp(int(y[i] / HC), 0, HH - 1);
			near.clear();
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
						if (j == i) continue;
						int bj = body[j];
						if (bj >= 0)
						{
							if (bj == body[i]) continue;
							if (std::find(near.begin(), near.end(), bj) == near.end()) near.push_back(bj);
							continue;
						}
						// j is a loose grain: grains collide as circles; grain vs body is handled from the grain's side
						if (body[i] != FREE || bj != FREE || j < i) continue;
						float pix = x[i], piy = y[i], pjx = x[j], pjy = y[j];
						float dx = pix - pjx, dy = piy - pjy, d2 = dx * dx + dy * dy;
						if (d2 >= D * D) continue;
						float d = std::sqrt(d2);
						if (d < 1e-6f) { dx = 0; dy = -1; d = 1; }
						float nx = dx / d, ny = dy / d, pen = D - d;
						float wi = 1.0f / dens[i], wj = 1.0f / dens[j], lam = pen / (wi + wj);
						Apply(i, pix, piy, nx * lam, ny * lam);
						Apply(j, pjx, pjy, -nx * lam, -ny * lam);
						if (recording) conts.push_back({ i, -1, j, nx, ny, 0, 0 });
					}
				}
			}
			for (int b : near) CollideWithBody(i, b);
		}
	}

	bool Solid(Env &env, int px, int py)
	{
		if (inGrid(px, py) && pendPix[size_t(py) * W + px]) return true;
		return env.Blocked(px, py);
	}
	void StaticCollide(Env &envRef, int i)
	{
		struct Wrap { World *w; Env &e; bool Blocked(int a, int b) { return w->Solid(e, a, b); } } env{ this, envRef };
		float px, py; Pos(i, px, py);
		int cx = int(std::floor(px + 0.5f)), cy = int(std::floor(py + 0.5f));
		if (env.Blocked(cx, cy))
		{
			// centre is inside an obstacle pixel: back out toward the side it came from
			const int dirs[4][2] = { {0,-1},{0,1},{-1,0},{1,0} };
			float best = 1e9f; int bx = 0, by = -1;
			for (auto &d : dirs)
			{
				if (env.Blocked(cx + d[0], cy + d[1])) continue;
				float tx = cx + d[0] - ox[i], ty = cy + d[1] - oy[i], dd = tx * tx + ty * ty;
				if (dd < best) { best = dd; bx = d[0]; by = d[1]; }
			}
			float tx, ty;
			if (best < 1e9f)
			{
				tx = bx ? cx + bx * (0.5f + R + 0.001f) : px;
				ty = by ? cy + by * (0.5f + R + 0.001f) : py;
			}
			else { tx = ox[i]; ty = oy[i]; }
			float dx = tx - px, dy = ty - py, dl = std::sqrt(dx * dx + dy * dy);
			if (dl > 1e-6f)
			{
				float nx = dx / dl, ny = dy / dl, w = InvMass(i, px, py, nx, ny);
				if (w > 0) { Apply(i, px, py, nx * dl / w, ny * dl / w); corr[i] = std::max(corr[i], dl); if (recording) conts.push_back({ i, -1, -1, nx, ny, 0, 0 }); }
			}
			return;
		}
		for (int yy = cy - 1; yy <= cy + 1; yy++)
			for (int xx = cx - 1; xx <= cx + 1; xx++)
			{
				if ((xx == cx && yy == cy) || !env.Blocked(xx, yy)) continue;
				float qx = std::clamp(px, xx - 0.5f, xx + 0.5f), qy = std::clamp(py, yy - 0.5f, yy + 0.5f);
				float dx = px - qx, dy = py - qy, d2 = dx * dx + dy * dy;
				if (d2 >= R * R || d2 < 1e-12f) continue;
				if (qx != px && qy != py)
				{
					int sx = px > qx ? 1 : -1, sy = py > qy ? 1 : -1;
					if (env.Blocked(xx + sx, yy) || env.Blocked(xx, yy + sy)) continue; // inner corner of a flat surface
				}
				float d = std::sqrt(d2), nx = dx / d, ny = dy / d, pen = R - d;
				float w = InvMass(i, px, py, nx, ny);
				if (w <= 0) return;
				float lam = pen / w;
				Apply(i, px, py, nx * lam, ny * lam);
				corr[i] = std::max(corr[i], pen);
				if (recording) conts.push_back({ i, -1, -1, nx, ny, 0, 0 });
				// friction against the ground
				float rxv = px - ox[i], ryv = py - oy[i], rn = rxv * nx + ryv * ny;
				float tx = rxv - rn * nx, ty = ryv - rn * ny, tl = std::sqrt(tx * tx + ty * ty);
				if (tl < 1e-7f) continue;
				tx /= tl; ty /= tl;
				float wt = InvMass(i, px, py, tx, ty);
				if (wt <= 0) continue;
				float lt = tl / wt;
				if (tl > fr[i] * pen * 1.5f) lt = std::min(lt, fr[i] * lam);
				Apply(i, px, py, -tx * lt, -ty * lt);
			}
	}

	// ---------------- fracture ----------------
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
		Split(b, &key, &dust, kd);
	}

	void Fractures()
	{
		int budget = 4;
		size_t nb = bodies.size();
		for (size_t b = 0; b < nb && budget > 0; b++)
		{
			Body &B = bodies[b];
			if (!B.alive || B.stat || B.asleep) continue;
			B.recent = B.recent * 0.5f + B.impact; B.impact = 0;
			if (B.cool > 0) { B.cool--; continue; }
			float weakest = 1e9f, od = -1; int origin = -1;
			for (int i : B.m)
			{
				weakest = std::min(weakest, str[i]);
				if (boundary[i] && corr[i] > od) { od = corr[i]; origin = i; }
			}
			float ratio = (B.recent / weakest) * (B.recent / weakest);
			if (ratio > 1 && origin >= 0) { budget--; B.recent = 0; Crack(int(b), origin, ratio); }
		}
		for (auto &B : bodies) B.impact = 0;
		for (int i : tracked) corr[i] = 0;
	}
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
		// strength = impact speed change (px/frame) the material survives
		set(PT_IRON, 7.8f, 0.45f, 20.0f, K_METAL);
		set(PT_METL, 7.8f, 0.45f, 19.0f, K_METAL);
		set(PT_GOLD, 19.3f, 0.40f, 16.0f, K_METAL);
		set(PT_WOOD, 0.6f, 0.65f, 9.0f, K_WOOD);
		set(PT_GLAS, 2.5f, 0.25f, 4.0f, K_GLASS);
		set(PT_BRCK, 2.0f, 0.75f, 5.5f, K_STONE);
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
		// velocity given by TPT (e.g. from other elements) becomes an impulse; TPT itself must not move rigid parts
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
		if (oxp >= 0 && oyp >= 0 && oxp < XRES && oyp < YRES && sim->pmap[oyp][oxp] && ID(sim->pmap[oyp][oxp]) == i)
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

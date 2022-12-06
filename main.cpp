#include <iostream>
#include <cstdint>
#include <cmath>

#include <string>
#include <ostream>

#include <algorithm>
#include <thread>
#include <mutex>

#include "BWindow/GDIWindow.h"


/*
var scene = 
	{
		gravity : -9.81,
		dt : 1.0 / 120.0,
		numIters : 100,
		frameNr : 0,
		overRelaxation : 1.9,
		obstacleX : 0.0,
		obstacleY : 0.0,
		obstacleRadius: 0.15,
		paused: false,
		sceneNr: 0,
		showObstacle: false,
		showStreamlines: false,
		showVelocities: false,	
		showPressure: false,
		showSmoke: true,
		fluid: null
	};*/

struct vec2 {
	float x, y;
};


template<size_t NUM_CELLS_X, size_t NUM_CELLS_Y>
struct FluidGrid {
	// float hor [NUM_CELLS_Y][NUM_CELLS_X + 1]{}; // horizontal velocities
	// float vert[NUM_CELLS_Y + 1][NUM_CELLS_X]{}; // vertical velocities
	float hor [NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{}; // horizontal velocities
	float vert[NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{}; // vertical velocities
	// float pressure[NUM_CELLS_Y][NUM_CELLS_X]{};

	enum Field {
		VEL_X,
		VEL_Y
	};

	inline void addVel(const float vHor, const float vVert) {
		for(size_t y = 0; y < NUM_CELLS_Y; y++)
			for(size_t x = 0; x < NUM_CELLS_X + 1; x++)
				hor[y][x] += vHor;
		for(size_t y = 0; y < NUM_CELLS_Y + 1; y++)
			for(size_t x = 0; x < NUM_CELLS_X; x++)
				vert[y][x] += vVert;
	}

	inline float divergence(const size_t x, const size_t y) const { // calculate divergence (outflow at specified cell)
		return hor[y][x + 1] - hor[y][x] // horizontal divergence
			+ vert[y + 1][x] - vert[y][x]; // vertical divergence
	}

	inline void solveDivergence(const size_t num_iter = 100, const float overRelaxation = 1.9f) { // change velocities so divergence becomes 0
		for(size_t i = 0; i < num_iter; i++) {
			for(size_t y = 0; y < NUM_CELLS_Y; y++) {
				for(size_t x = 0; x < NUM_CELLS_X; x++) {
					const float currentDivergence = divergence(x, y);
					const float correctionFactor = -currentDivergence / 4.f * overRelaxation; // TODO: adapt for walls (replace 4)

					hor[y][x + 1]	-= correctionFactor;
					hor[y][x]		+= correctionFactor;
					vert[y + 1][x]	-= correctionFactor;
					vert[y][x]		+= correctionFactor;
				}
			}
		}
	}

	inline float sample(const Field field, float x, float y) const {
		// x -= .5f;
		// y -= .5f;
		// if(x < 0) x = 0;
		// if(y < 0) y = 0;

		// const size_t x0 = (size_t)x;
		// const size_t y0 = (size_t)y;
		// const float xt = x - x0;
		// const float yt = y - y0;

		// const float vHor0 = hor[y0    ][x0] * (1.f - xt) + hor[y0    ][x0 + 1] * xt;
		// const float vHor1 = hor[y0 + 1][x0] * (1.f - xt) + hor[y0 + 1][x0 + 1] * xt;

		// const float vVert0 = vert[y0][x0    ] * (1.f - yt) + vert[y0 + 1][x0    ] * yt;
		// const float vVert1 = vert[y0][x0 + 1] * (1.f - yt) + vert[y0 + 1][x0 + 1] * yt;

		// return {
		// 	vHor0 * (1.f - yt) + vHor1 * yt,
		// 	vVert0 * (1.f - xt) + vVert1 * xt
		// };

		// var n = this.numY;
		// var h = this.h;
		const float h = 1.f;
		const float h1 = 1.0 / h;
		const float h2 = 0.5 * h;

		x = std::max<float>(std::min<float>(x, NUM_CELLS_X), 1);
		y = std::max<float>(std::min<float>(y, NUM_CELLS_Y), 1);

		float dx = 0.f;
		float dy = 0.f;

		// const float (*f)[NUM_CELLS_Y + 1][NUM_CELLS_X + 1];
		const float (*f)[NUM_CELLS_X + 1];

		switch (field) {
			// case U_FIELD: f = this.u; dy = h2; break;
			// case V_FIELD: f = this.v; dx = h2; break;
			// case S_FIELD: f = this.m; dx = h2; dy = h2; break
			case VEL_X: f = hor; dy = h2; break;
			case VEL_Y: f = vert; dx = h2; break;
		}

		const size_t x0 = std::min<size_t>(std::floor((x-dx)*h1), NUM_CELLS_X-1);
		const float tx = ((x-dx) - x0*h) * h1;
		const size_t x1 = std::min<size_t>(x0 + 1, NUM_CELLS_X-1);
		
		const size_t y0 = std::min<size_t>(std::floor((y-dy)*h1), NUM_CELLS_Y-1);
		const float ty = ((y-dy) - y0*h) * h1;
		const size_t y1 = std::min<size_t>(y0 + 1, NUM_CELLS_Y-1);

		const float sx = 1.0 - tx;
		const float sy = 1.0 - ty;

		const float val =
			sx*sy * f[y0][x0] +
			tx*sy * f[y0][x1] +
			tx*ty * f[y1][x1] +
			sx*ty * f[y1][x0];
		// const float val = sx*sy * f[x0*n + y0] +
		// 	tx*sy * f[x1*n + y0] +
		// 	tx*ty * f[x1*n + y1] +
		// 	sx*ty * f[x0*n + y1];
		
		return val;
	}

	inline float avgVelX(const size_t x, const size_t y) const {
		return (hor[y][x] + hor[y][x + 1]
			+ hor[y + 1][x] + hor[y + 1][x + 1]) * .25f;
	}

	inline float avgVelY(const size_t x, const size_t y) const {
		return (vert[y][x] + vert[y + 1][x]
			+ vert[y][x + 1] + vert[y + 1][x + 1]) * .25f;
	}

	static inline void update(FluidGrid& prev, FluidGrid& next, const float dt) {
		memcpy(next.hor, prev.hor, sizeof(float) * (NUM_CELLS_X + 1) * NUM_CELLS_Y);
		memcpy(next.vert, prev.vert, sizeof(float) * NUM_CELLS_X * (NUM_CELLS_Y + 1));

		for(size_t y = 1; y < NUM_CELLS_Y; y++) {
			for(size_t x = 1; x < NUM_CELLS_X; x++) {
				// update horizontal velocities:
				{
					float px = x;
					float py = y + .5f;
					const float velX = prev.hor[y][x];
					const float velY = prev.avgVelY(x, y);
					px -= dt*velX;
					py -= dt*velY;
					const float newVelX = prev.sample(VEL_X, px, py);
					next.hor[y][x] = newVelX;
				}

				// update vertical velocities:
				{
					float px = x + .5f;
					float py = y;
					const float velX = prev.avgVelX(x, y);
					const float velY = prev.vert[y][x];
					px -= dt*velX;
					py -= dt*velY;
					const float newVelY = prev.sample(VEL_Y, px, py);
					next.vert[y][x] = newVelY;
				}
			}
		}

		prev = next;
	}
};

int main() {
	GDIWindow win(800, 800);

	static constexpr size_t CELLS_X = 100;
	static constexpr size_t CELLS_Y = 100;

	// TODO: no memory leak :)
	FluidGrid<CELLS_X, CELLS_Y>* vCurrent = new FluidGrid<CELLS_X, CELLS_Y>;
	FluidGrid<CELLS_X, CELLS_Y>* vNext = new FluidGrid<CELLS_X, CELLS_Y>;


	// float density = 1000.f;

	while(!win.shouldClose()) {
		vCurrent->addVel(0, -9.81f);
		FluidGrid<CELLS_X, CELLS_Y>::update(*vCurrent, *vNext, .016f);

		for(size_t y = 0; y < win.height; y++) {
			for(size_t x = 0; x < win.width; x++) {
				// win.graphics.setPixel(x, y, 0xFF00FF00);
				win.graphics.setPixel(x, y, (x % 100 < 50) ? 0xFF00FF00 : 0xFFFF0000);
			}
		}

		win.updateScreen();

		win.pollMsg();
	}

	return 0;
}
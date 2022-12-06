#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
// #include <thread>
// #include <mutex>

#include "BWindow/GDIWindow.h"

#include "windows.h"


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


enum class Field : uint8_t {
	VEL_X,
	VEL_Y,
	SMOKE
};

template<size_t NUM_CELLS_X, size_t NUM_CELLS_Y>
struct FluidGrid {
	float hor  [NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{}; // horizontal velocities
	float vert [NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{}; // vertical velocities
	float s    [NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{}; // scaling (0 = wall, 1 = fluid)
	float smoke[NUM_CELLS_Y + 1][NUM_CELLS_X + 1]{};
	// float pressure[NUM_CELLS_Y][NUM_CELLS_X]{};

	inline FluidGrid() {
		for(size_t y = 0; y < NUM_CELLS_Y + 1; y++)
			for(size_t x = 0; x < NUM_CELLS_X + 1; x++)
			// for(size_t x = 1; x < NUM_CELLS_X; x++)
				s[y][x] = 1.f;
	}

	inline void addVel(const float vHor, const float vVert) {
		for(size_t y = 1; y < NUM_CELLS_Y; y++)
			for(size_t x = 1; x < NUM_CELLS_X + 1; x++)
				hor[y][x] += vHor;
		for(size_t y = 1; y < NUM_CELLS_Y + 1; y++)
			for(size_t x = 1; x < NUM_CELLS_X; x++)
				vert[y][x] += vVert;
	}

	inline float divergence(const size_t x, const size_t y) const { // calculate divergence (outflow at specified cell)
		return hor[y][x + 1] - hor[y][x] // horizontal divergence
			+ vert[y + 1][x] - vert[y][x]; // vertical divergence
	}

	// inline void solveDivergence(const size_t num_iter = 100, const float overRelaxation = 1.9f) { // change velocities so divergence becomes 0
	// 	for(size_t i = 0; i < num_iter; i++) {
	// 		for(size_t y = 0; y < NUM_CELLS_Y; y++) {
	// 			for(size_t x = 0; x < NUM_CELLS_X; x++) {
	// 				const float currentDivergence = divergence(x, y);
	// 				const float correctionFactor = -currentDivergence / 4.f * overRelaxation; // TODO: adapt for walls (replace 4)

	// 				hor[y][x + 1]	-= correctionFactor;
	// 				hor[y][x]		+= correctionFactor;
	// 				vert[y + 1][x]	-= correctionFactor;
	// 				vert[y][x]		+= correctionFactor;
	// 			}
	// 		}
	// 	}
	// }

	inline void solveDivergence(const size_t num_iter = 100, const float overRelaxation = 1.9f) { // change velocities so divergence becomes 0
		for(size_t i = 0; i < num_iter; i++) {
			for(size_t y = 1; y < NUM_CELLS_Y; y++) {
				for(size_t x = 1; x < NUM_CELLS_X; x++) {
					if(s[y][x] == 0.f) continue; // inside wall

					const float sx0 = s[y][x-1];
					const float sx1 = s[y][x+1];
					const float sy0 = s[y-1][x];
					const float sy1 = s[y+1][x];
					const float s = sx0 + sx1 + sy0 + sy1;

					if (s == 0.f) continue; // surrounded by walls

					const float currentDivergence = divergence(x, y);
					const float correctionFactor = -currentDivergence / s * overRelaxation; // TODO: adapt for walls (replace 4)

					hor[y][x  ] -= correctionFactor * sx0;
					hor[y][x+1] += correctionFactor * sx1;
					vert[y  ][x] -= correctionFactor * sy0;
					vert[y+1][x] += correctionFactor * sy1;
					// hor[y][x + 1]	-= correctionFactor;
					// hor[y][x]		+= correctionFactor;
					// vert[y + 1][x]	-= correctionFactor;
					// vert[y][x]		+= correctionFactor;
				}
			}
		}
	}

	inline void extrapolate() {
		for (size_t x = 0; x < NUM_CELLS_X; x++) {
			hor[0][x] = hor[1][x];
			hor[NUM_CELLS_Y][x] = hor[NUM_CELLS_Y - 1][x]; 
		}
		for (size_t y = 0; y < NUM_CELLS_Y; y++) {
			vert[y][0] = vert[y][1];
			vert[y][NUM_CELLS_X] = vert[y][NUM_CELLS_X - 1];
		}
	}

	inline float sample(const Field field, float x, float y) const {
		const float h = 1.f;
		const float h1 = 1.f / h;
		const float h2 = .5f * h;

		x = std::max<float>(std::min<float>(x, NUM_CELLS_X + 1), 1);
		y = std::max<float>(std::min<float>(y, NUM_CELLS_Y + 1), 1);

		float dx = 0.f;
		float dy = 0.f;

		// const float (*f)[NUM_CELLS_Y + 1][NUM_CELLS_X + 1];
		const float (*f)[NUM_CELLS_X + 1];

		switch (field) {
			// case U_FIELD: f = this.u; dy = h2; break;
			// case V_FIELD: f = this.v; dx = h2; break;
			// case S_FIELD: f = this.m; dx = h2; dy = h2; break
			case Field::VEL_X: f = hor; dy = h2; break;
			case Field::VEL_Y: f = vert; dx = h2; break;
			case Field::SMOKE: f = smoke; dx = h2; dy = h2; break;
			// case S_FIELD: f = this.m; dx = h2; dy = h2; break
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
		// memcpy(next.hor, prev.hor, sizeof(float) * (NUM_CELLS_X + 1) * NUM_CELLS_Y);
		// memcpy(next.vert, prev.vert, sizeof(float) * NUM_CELLS_X * (NUM_CELLS_Y + 1));

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
					const float newVelX = prev.sample(Field::VEL_X, px, py);
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
					const float newVelY = prev.sample(Field::VEL_Y, px, py);
					next.vert[y][x] = newVelY;
				}
			}
		}

		// prev = next;
		memcpy(prev.hor, next.hor, sizeof(float) * (NUM_CELLS_X + 1) * (NUM_CELLS_Y + 1));
		memcpy(prev.vert, next.vert, sizeof(float) * (NUM_CELLS_X + 1) * (NUM_CELLS_Y + 1));
	}

	static inline void updateSmoke(FluidGrid& prev, FluidGrid& next, const float dt) {
		for(size_t y = 1; y < NUM_CELLS_Y; y++) {
			for(size_t x = 1; x < NUM_CELLS_X; x++) {
				if(prev.s[y][x] == 0) continue;

				const float velX = (prev.hor[y][x] + prev.hor[y][x + 1]) * .5f;
				const float velY = (prev.vert[y][x] + prev.vert[y + 1][x]) * .5f;

				// const float sourceX = x + .5f - velX * dt;
				// const float sourceY = y + .5f - velY * dt;

				const float sourceX = x + .50f - velX * dt;
				const float sourceY = y + .50f - velY * dt;

				next.smoke[y][x] = prev.sample(Field::SMOKE, sourceX, sourceY);
			}
		}

		memcpy(prev.smoke, next.smoke, sizeof(float) * (NUM_CELLS_X + 1) * (NUM_CELLS_Y + 1));
	}
};

int main() {
	GDIWindow win(800, 800);

	const auto drawDir = [&win](const vec2& pos, const vec2& dir, const float size){
			win.graphics.line((int)pos.x, (int)pos.y, (int)(pos.x + dir.x * size), (int)(pos.y + dir.y * size), 0xFFFFFFFF);
			win.graphics.setPixel((int)pos.x, (int)pos.y, 0xFF0000FF);
		};

	static constexpr size_t CELLS_X = 400;
	static constexpr size_t CELLS_Y = 400;

	// TODO: no memory leak :)
	FluidGrid<CELLS_X, CELLS_Y>* vCurrent = new FluidGrid<CELLS_X, CELLS_Y>;
	FluidGrid<CELLS_X, CELLS_Y>* vNext = new FluidGrid<CELLS_X, CELLS_Y>;

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq); // TODO


	LARGE_INTEGER lastTime;
	QueryPerformanceCounter(&lastTime);


	for(size_t x = 0; x < CELLS_X + 1; x++) {
		vCurrent->s[0][x] = 0;
		vCurrent->s[CELLS_Y][x] = 0;
	}
	for(size_t y = 0; y < CELLS_Y + 1; y++)
		vCurrent->s[y][0] = 0;

	for(size_t y = 0; y < CELLS_Y / 8; y++) {
		for(size_t x = 0; x < CELLS_X / 8; x++) {
			vCurrent->s[CELLS_Y / 3 * 1 + y][CELLS_X / 5 + x] = 0;
			vCurrent->s[CELLS_Y / 3 * 2 - y][CELLS_X / 5 + x] = 0;
		}
	}

	// float density = 1000.f;
	while(!win.shouldClose()) {
		LARGE_INTEGER currentTime;
		QueryPerformanceCounter(&currentTime);

		const double dt = (currentTime.QuadPart - lastTime.QuadPart) * 1. / freq.QuadPart;
		lastTime = currentTime;

		std::cout << "dt: " << dt << "s\n";

		const size_t NUM_SMOKE_TRAILS = 10;
		for(size_t y = CELLS_Y / NUM_SMOKE_TRAILS / 2; y < CELLS_Y; y += CELLS_Y / NUM_SMOKE_TRAILS)
			vCurrent->smoke[y][0] = 3.f;

		for(size_t y = 1; y < CELLS_Y - 3; y++)
			vCurrent->hor[y][1] = 30;
		// for(size_t y = 1; y < CELLS_Y - 2; y++)
		// 	vCurrent->vert[y][1] = 0;

		// vCurrent->hor[54][5] = 50;

		// vCurrent->addVel(0, 9.81f);

		// static bool first=true;
		// if(first ||GetAsyncKeyState('S') & 0x8000)
		vCurrent->solveDivergence(10, 1.9f);
		// first=false;

		// const size_t iter = 1;
		// for(size_t i = 0; i < iter; i++)
		// 	FluidGrid<CELLS_X, CELLS_Y>::update(*vCurrent, *vNext, dt / iter);

		vCurrent->extrapolate();

		FluidGrid<CELLS_X, CELLS_Y>::update(*vCurrent, *vNext, dt);

		FluidGrid<CELLS_X, CELLS_Y>::updateSmoke(*vCurrent, *vNext, dt);


		// --- graphics:
		for(size_t y = 0; y < win.height; y++) {
			for(size_t x = 0; x < win.width; x++) {
				const size_t cellX = x * (CELLS_X + 1) / win.width;
				const size_t cellY = y * (CELLS_Y + 1) / win.height;

				static constexpr auto clamp = [](float v, float min, float max){
						return (v < min) ? min : (v > max) ? max : v;
					};

				const float wall = vCurrent->s[cellY][cellX];

				const float smoke = 1.f - clamp(vCurrent->sample(Field::SMOKE, cellX, cellY), 0.f, 1.f);

				const uint8_t vX = smoke * wall * (clamp(vCurrent->sample(Field::VEL_X, cellX, cellY) * .01f, -.5f, .5f) * 255 + 128);
				const uint8_t vY = smoke * wall * (clamp(vCurrent->sample(Field::VEL_Y, cellX, cellY) * .01f, -.5f, .5f) * 255 + 128);

				win.graphics.setPixel(x, y,
					0xFF << 24
					| vX << 16
					| vY << 8
				);
				// win.graphics.setPixel(x, y, (x % 100 < 50) ? 0xFF00FF00 : 0xFFFF0000);
			}
		}


		if(GetAsyncKeyState('L') & 0x8000) {
			for(size_t y = 0; y < CELLS_Y; y++) {
				for(size_t x = 0; x < CELLS_X; x++) {
					drawDir(
						vec2 {
							(x + .5f) * win.width / CELLS_X,
							(y + .5f) * win.height / CELLS_Y
						},
						vec2 {
							vCurrent->sample(Field::VEL_X, x, y),
							vCurrent->sample(Field::VEL_Y, x, y)
						},
						5.f);
				}
			}
	
			// for(size_t y = 0; y < CELLS_Y; y++) {
			// 	for(size_t x = 0; x < CELLS_X + 1; x++) {
			// 		const int startX = x * win.width / CELLS_X;
			// 		const int startY = (y + .5f) * win.height / CELLS_Y;
			// 		win.graphics.line(
			// 			startX,
			// 			startY,
			// 			startX + std::min<float>(vCurrent->hor[y][x] * 5.f, 15),
			// 			startY,
			// 			0xFFFFFFFF
			// 		);
			// 		win.graphics.fillCircle(startX, startY, 1, 0xFFFF00FF);
			// 	}
			// }
			// for(size_t y = 0; y < CELLS_Y + 1; y++) {
			// 	for(size_t x = 0; x < CELLS_X; x++) {
			// 		const int startX = (x + .5f) * win.width / CELLS_X;
			// 		const int startY = y * win.height / CELLS_Y;
			// 		win.graphics.line(
			// 			startX,
			// 			startY,
			// 			startX,
			// 			startY + std::min<float>(vCurrent->vert[y][x] * 5.f, 15),
			// 			0xFFFFFFFF
			// 		);
			// 		win.graphics.fillCircle(startX, startY, 1, 0xFFFF00FF);
			// 	}
			// }
		}

		win.updateScreen();

		win.pollMsg();
	}

	return 0;
}
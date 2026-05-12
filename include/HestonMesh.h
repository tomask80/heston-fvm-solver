#pragma once
class HestonMesh
{
public:
	int Nx, Ny;
	double hx, hy;
	double X_l, y_min = 0.0;

	HestonMesh(int nx, int ny, double xl, double xr, double ymax) // asi zabdunuty xr
		: Nx(nx), Ny(ny), X_l(xl) {
		hx = (xr - xl) / Nx;
		hy = ymax / Ny;
	}

	inline int getIndex(int i, int j) const {
		return j * Nx + i;
	}

	double getX(int i) const { return X_l + (i + 0.5) * hx; } // v sulade s clankom 
	double getY(int j) const { return (j + 0.5) * hy; }
};


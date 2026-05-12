#pragma once
#include "HestonMesh.h"
#include "iostream"
#include "Eigen/sparseLu"
#include "HestonParams.h"
#include "HestonAnalytical.h"
#include <fstream>

class HestonSolver
{
private:
	HestonParams params;
	HestonMesh mesh;
	Eigen::SparseMatrix<double> L; // matica sustavy
	Eigen::VectorXd u; // akutalne riesenie (cena opcie)
	std::vector<double> dirichlet_right_coeffs; // pre pravu stranu dir OP lebo je zavisla od casu

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	double current_tau = 0.0;

	// Komponenty tenzora difuzie B podla (15)
	// B = 0.5 * y * [ {1, rho*sigma }, {rho*sigma, sigma*sigma} ]
	double getB11(double y) { return 0.5 * y; }
	double getB12(double y, double rho, double sigma) { return 0.5 * y * rho * sigma; }
	double getB22(double y, double sigma) { return 0.5 * y * sigma * sigma; }

	// Komponenty advekcneho vektora A podla (16)
	double getA1(double y, double r, double rho, double sigma) {
		return -(r - 0.5 * y - 0.5 * rho * sigma);
	}
	double getA2(double y, double kappa, double theta, double sigma, double lambda = 0.0) {
		return -(kappa * (theta - y) - lambda * y - 0.5 * sigma * sigma);
	}

public:
	HestonSolver(HestonParams p, int nx, int ny) : params(p), mesh(nx, ny, p.X_l, p.X_r, p.y_maxl) {
		u.resize(mesh.Nx * mesh.Ny);
		L.resize(mesh.Nx * mesh.Ny, mesh.Nx * mesh.Ny);
	}

	// Accessors for mesh spacing so caller can set dt = hx * hy
	double getHx() const { return mesh.hx; }
	double getHy() const { return mesh.hy; }

	void setInitialCondition();
	double calculateCpq(int i, int j, int di, int dj,double scale_hx = 1.0);
	void buildMatrix(double dt);
	void step(double dt); // jeden casovy krok (implicit Euler) 
	void initializeSolver();

	// optional: export helpers (used in main)
	void exportToCSV(const std::string& filename);
	void exportToRealPrices(const std::string& filename, double E);

	void exportFromAnalytical(const std::string& filename,HestonAnalytical& exact);

	// Add accessors for u and mesh
	const Eigen::VectorXd& getU() const { return u; }
	const HestonMesh& getMesh() const { return mesh; }
};


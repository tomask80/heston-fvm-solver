#include "HestonSolver.h"

void HestonSolver::setInitialCondition()
{
	for (size_t j = 0; j < mesh.Ny; j++)
	{
		for (size_t i = 0; i < mesh.Nx; i++)
		{
			double x = mesh.getX(i);
			int idx = mesh.getIndex(i,j);
			// vypocet transformovaneho payoffu: max(e^x -1,0)
			u[idx] = std::max(std::exp(x) - 1.0, 0.0);
		}
	}
	std::cout << "Initial condition set." << std::endl;
}

double HestonSolver::calculateCpq(int i, int j, int di, int dj,double scale_hx) {
    double hx = mesh.hx;
    double hy = mesh.hy;

    // suradnice y v stredoch bunky a prislusnycg hran
    double y_p = mesh.getY(j);               // Pre hrany pe, pw
    double y_n = mesh.getY(j) + 0.5 * hy;    // Pre hranu pn
    double y_s = mesh.getY(j) - 0.5 * hy;    // Pre hranu ps

    // --- 1. vycislenie tenzora B v stredoch hran ---
    double b11_pe = getB11(y_p);
    double b11_pw = getB11(y_p); 

    double b22_pn = getB22(y_n, params.sigma);
    double b22_ps = getB22(y_s, params.sigma);

    double b12_pe = getB12(y_p, params.rho, params.sigma);
    double b12_pw = getB12(y_p, params.rho, params.sigma);
	double b21_pn = getB12(y_n, params.rho, params.sigma); // b21 = b12 symetrra tenzora    
    double b21_ps = getB12(y_s, params.rho, params.sigma);

    // --- 2. vycislenie vektora A v stredoch hran ---
    double a1_pe = getA1(y_p, params.r, params.rho, params.sigma);
    double a1_pw = getA1(y_p, params.r, params.rho, params.sigma);
    double a2_pn = getA2(y_n, params.kappa, params.theta, params.sigma);
    double a2_ps = getA2(y_s, params.kappa, params.theta, params.sigma);

    // --- 3. vypocet koeficientov a_pq (39) ---
    double a_pe = -0.5 * hy * a1_pe;
    double a_pw = 0.5 * hy * a1_pw;
    double a_pn = -0.5 * hx * a2_pn;
    double a_ps = 0.5 * hx * a2_ps;

    // --- 4. difuzne koeficienty (40 a 41) ---
    // scale_hx upravuje efektivnu vzdialenost pre x-gradient pri hranicnych bunkach.
    // Pre vnutorne bunky scale_hx = 1.0  =>  hx_eff = hx   (vzdialenost stred-stred)
    // Pre hranicne bunky scale_hx = 0.5  =>  hx_eff = hx/2 (vzdialenost stred-hrana)
    // Pomer hy/hx_eff = hy/(scale_hx * hx), cize delime scale_hx.
    double b_pe = (hy / (scale_hx * hx)) * b11_pe + (b21_pn / 4.0) - (b21_ps / 4.0);
    double b_pw = (hy / (scale_hx * hx)) * b11_pw - (b21_pn / 4.0) + (b21_ps / 4.0);
    double b_pn = (hx / hy) * b22_pn + (b12_pe / 4.0) - (b12_pw / 4.0);
    double b_ps = (hx / hy) * b22_ps - (b12_pe / 4.0) + (b12_pw / 4.0);

    // --- 5. vypocet diagonalnych koeficientov (42) ---
    double b_pne = (b12_pe / 4.0) + (b21_pn / 4.0);
    double b_psw = (b12_pw / 4.0) + (b21_ps / 4.0);
    double b_pnw = -(b12_pw / 4.0) - (b21_pn / 4.0);
    double b_pse = -(b12_pe / 4.0) - (b21_ps / 4.0);

    // --- 6. navrat hodnoty c_pq podla smeru (di, dj) ---
    // Pre priamych susedov vraciame a_pq + b_pq
    if (di == 1 && dj == 0) return a_pe + b_pe; // East
    else if (di == -1 && dj == 0) return a_pw + b_pw; // West
    else if (di == 0 && dj == 1) return a_pn + b_pn; // North
    else if (di == 0 && dj == -1) return a_ps + b_ps; // South

    // Pre diagonalnych susedov vraciame len b_pq (a_pq je tam nulove)
    else if (di == 1 && dj == 1) return b_pne;       // North-East
    else if (di == -1 && dj == 1) return b_pnw;       // North-West
    else if (di == 1 && dj == -1) return b_pse;       // South-East
    else if (di == -1 && dj == -1) return b_psw;       // South-West

    return 0.0;
}

void HestonSolver::buildMatrix(double dt) {
    std::vector<Eigen::Triplet<double>> triplets;
    double dt_over_mp = dt / (mesh.hx * mesh.hy);

    // pripravime si vektor pre koeficienty pravej Dirichletovej podmienky
    dirichlet_right_coeffs.assign(mesh.Nx * mesh.Ny, 0.0);

    for (int j = 0; j < mesh.Ny; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            int p = mesh.getIndex(i, j);
            double sum_cpq = 0.0;

            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    if (di == 0 && dj == 0) continue;

                    int ni = i + di;
                    int nj = j + dj;

                    // --- ZRKADLENIE (Ghost Cells pre y-hranice) ---
                    int nj_eff = nj;
                    if (nj < 0) nj_eff = 0;
                    if (nj >= mesh.Ny) nj_eff = mesh.Ny - 1;

                    // --- DIRICHLET na x-hraniciach ---
                    if (ni < 0) {
                        // Lava hranica: u_W = 0 predpisane na hrane x = X_l.
                        // Vzdialenost od stredu bunky i=0 po hranu je hx/2, preto scale_hx = 0.5.
                        // Platí len pre priameho zapadneho suseda (di=-1, dj=0).
                        // Diagonalni susedia (di=-1, dj=+-1) maju x-gradient cez hx/2 tiez,
                        // ale ich b_pnw/b_psw nezavisia od hx (len od b12, b21) takze scale nema efekt.
                        double cpq = calculateCpq(i, j, di, dj, 0.5);
                        sum_cpq += cpq;
                    }
                    else if (ni >= mesh.Nx) {
                        // Prava hranica: u_E = bc_right_val predpisane na hrane x = X_r.
                        // Vzdialenost od stredu bunky i=Nx-1 po hranu je hx/2, preto scale_hx = 0.5.
                        double cpq = calculateCpq(i, j, di, dj, 0.5);
                        sum_cpq += cpq;
                        dirichlet_right_coeffs[p] += dt_over_mp * cpq;
                    }
                    else {
                        // Vnutorna bunka — standardna vzdialenost hx medzi stredmi, scale_hx = 1.0
                        double cpq = calculateCpq(i, j, di, dj, 1.0);
                        int q_eff = mesh.getIndex(ni, nj_eff);
                        sum_cpq += cpq;
                        triplets.push_back(Eigen::Triplet<double>(p, q_eff, -dt_over_mp * cpq));
                    }
                }
            }

            // Rovnica 38: Pridanie finálneho diagonálneho prvku
            double diag = 1.0 + dt * params.r + dt_over_mp * sum_cpq;
            triplets.push_back(Eigen::Triplet<double>(p, p, diag));
        }
    }

    // Zostavenie matice - Eigen automaticky sčíta duplikáty (napr. zrkadlené bunky)
    L.setFromTriplets(triplets.begin(), triplets.end());
    L.makeCompressed();
 /*   for (int i = 0; i < mesh.Ny; i++)
    {
        for (size_t j = 0; j < mesh.Nx; j++)
        {
            std::cout << L.coeff(i, j) << " ";
        }
        std::cout << std::endl;
    }*/
	

    std::cout << "Matica L uspesne zostavena podla Diamond-cell schemy." << std::endl;
}


void HestonSolver::initializeSolver() {
    // Analýza a faktorizácia matice. 
    // Toto urobíme len RAZ pred spustením časovej slučky, čo obrovsky šetrí výkon.
    solver.analyzePattern(L);
    solver.factorize(L);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Chyba: Dekompozícia matice L zlyhala!" << std::endl;
        exit(1);
    }
    std::cout << "Matica L bola uspesne faktorizovana." << std::endl;
}

void HestonSolver::step(double dt) {
    current_tau += dt;

    // 1. Výpočet hodnoty Dirichletovej podmienky na pravej hranici (x = X_r)
    // u(X_r, y, tau) = exp(X_r) - exp(-r * tau)
    double bc_right_val = std::exp(params.X_r) - std::exp(-params.r * current_tau);

    // 2. Príprava vektora pravej strany (b). 
    // Na začiatku kroku je to len riešenie z predchádzajúceho času (u^{n-1}).
    Eigen::VectorXd b = u;

    // 3. Pridanie časovo závislých príspevkov z okrajov (využijeme vektor z buildMatrix)
    for (int idx = 0; idx < b.size(); ++idx) {
        b(idx) += dirichlet_right_coeffs[idx] * bc_right_val;
    }

    // 4. Riešenie systému: L * u^n = b
    u = solver.solve(b);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Chyba: Riesenie linearneho systemu zlyhalo v case tau = " << current_tau << std::endl;
    }
}


void HestonSolver::exportToCSV(const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Chyba: Nepodarilo sa otvorit subor pre zapis!" << std::endl;
        return;
    }

    for (int j = 0; j < mesh.Ny; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            double x = mesh.getX(i);
            double y = mesh.getY(j);
            double price = u[mesh.getIndex(i, j)];

            // Zápis v tvare: x, y, u
            file << x << "," << y << "," << price << "\n";
        }
    }
    file.close();
    std::cout << "Data uspesne exportovane do: " << filename << std::endl;
}

void HestonSolver::exportToRealPrices(const std::string& filename, double E) {
    std::ofstream file(filename);
    //file << "S,V,Price\n";

    for (int j = 0; j < mesh.Ny; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            double x = mesh.getX(i);
            double y = mesh.getY(j);
            double u_val = u[mesh.getIndex(i, j)];

            // Inverzná transformácia
            double S = E * std::exp(x);
            double realPrice = u_val * E;

            file << S << "," << y << "," << realPrice << "\n";
        }
    }
    file.close();
}

void HestonSolver::exportFromAnalytical(const std::string& filename, HestonAnalytical& exact)
{
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Chyba: Nepodarilo sa otvorit subor pre zapis!" << std::endl;
        return;
    }

    for (int j = 0; j < mesh.Ny; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            double x = mesh.getX(i);
            double y = mesh.getY(j);
            double price = exact.getCallPrice(std::exp(x), y, params.T, 1.0);

            // Zápis v tvare: x, y, u
            file << x << "," << y << "," << price << "\n";
        }
    }
    file.close();
    std::cout << "Data uspesne exportovane do: " << filename << std::endl;
}


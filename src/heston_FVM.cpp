#include <iostream>
#include <Eigen/Sparse>
#include "HestonSolver.h"
#include "HestonAnalytical.h"
#include <vector>
#include <cmath>
#include <math.h>

int main()
{
    HestonParams params;
    params.kappa = 5.0;
    params.theta = 0.07;
    params.sigma = 0.5;
    params.rho = -0.5;
    params.r = 0.1;
    params.T = 0.05;     // cas do splatnosti
    params.X_l = -7.0;   // log(S/E) min
    params.X_r = 3.0;    // log(S/E) max
    params.y_maxl = 1.0;  // Volatilita max

    int Nx = 40;
    int Ny = 20;

    HestonSolver solver(params, Nx, Ny);
    double hx = solver.getHx();
    double hy = solver.getHy();
    double dt = hx * hy;

    // Pre zachovanie pokrytia času T nastavíme počet krokov na ceil(T/dt)
    int N_ts = static_cast<int>(std::ceil(params.T / dt));
    if (N_ts <= 0) N_ts = 1;

    std::cout << "Using dt = hx*hy = " << dt << ", N_ts = " << N_ts << std::endl;

    solver.setInitialCondition();  // Nastaví payoff (tau = 0)
    solver.buildMatrix(dt);        // Vytvorí maticu L

    //presne riesenie
  /*  HestonAnalytical exact(params.kappa, params.theta, params.sigma, params.rho, params.r);
    double exactSol = exact.getCallPrice(exp(1.0), 0.5, params.T, 1.0);
    std::cout << "Exact sol = " << exactSol << std::endl;
    */
    
    solver.initializeSolver();     // Faktorizuje matic

 /*   std::cout << "Zacinam casovu slucku..." << std::endl;
    for (int n = 1; n <= N_ts; ++n) {
        solver.step(dt);
        std::cout << "Krok " << n << "/" << N_ts << " (tau = " << n * dt << ") dokonceny." << std::endl;
    }

    std::cout << "Vypocet uspesne ukonceny!" << std::endl;*/

    //L2 error
    //HestonAnalytical exact(params.kappa, params.theta, params.sigma, params.rho, params.r);
    //double total_l2_sum_sq = 0.0;

    //for (int n = 1; n <= N_ts; ++n) {
    //    solver.step(dt);
    //    double current_tau = n * dt; 

    //    double spatial_error_step = 0.0;
    //    for (int j = 0; j < Ny; ++j) {
    //        for (int i = 0; i < Nx; ++i) {
    //            double x = solver.getMesh().getX(i);
    //            double y = solver.getMesh().getY(j);

    //            double u_num = solver.getU()[solver.getMesh().getIndex(i, j)];
    //            double u_exact = exact.getCallPrice(std::exp(x), y, current_tau, 100.0); // tu som uz buchol E =100

    //            double diff = u_num - u_exact;
    //            spatial_error_step += (diff * diff) * hx * hy;
    //        }
    //    }
    //    total_l2_sum_sq += spatial_error_step * dt;

    //    if (n % 10 == 0) std::cout << "Step " << n << "/" << N_ts << " processed..." << std::endl;
    //}
    //double final_l2_error = std::sqrt(total_l2_sum_sq);
    //std::cout << "Celkova L2(I, Omega) chyba = " << final_l2_error << std::endl;

    HestonAnalytical exact(params.kappa, params.theta, params.sigma, params.rho, params.r);
    double total_l2_sum_sq = 0.0;
    double E_real = 1.0; 

    for (int n = 1; n <= N_ts; ++n) {
        solver.step(dt);
        double current_tau = n * dt;

        double spatial_error_step = 0.0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = solver.getMesh().getX(i);
                double y = solver.getMesh().getY(j);

                if (x >= -1.0 && x <= 1.0 && y >= 0.0 && y <= 1.0) {

                    double u_num_dimensionless = solver.getU()[solver.getMesh().getIndex(i, j)];
                    double u_num_real = u_num_dimensionless * E_real;

                    double S_real = E_real * std::exp(x);
                    double u_exact_real = exact.getCallPrice(S_real, y, current_tau, E_real);

                    double diff = u_num_real - u_exact_real;
                    spatial_error_step += (diff * diff) * hx * hy;
                }
            }
        }
        total_l2_sum_sq += spatial_error_step * dt;

        if (n % 10 == 0 || n == N_ts) {
            std::cout << "Step " << n << "/" << N_ts << " processed..." << std::endl;
        }
    }

    // 4. Finálna odmocnina
    double final_l2_error = std::sqrt(total_l2_sum_sq);
    std::cout << "Celkova L2(I, Omega_obs) chyba (E=100, x in [-1,1]) = " << final_l2_error << std::endl;


    //solver.exportToCSV("heston_solution.csv");
    //solver.exportToRealPrices("heston_real.csv", 100.0);
    //solver.exportFromAnalytical("heston_exact.csv", exact);

	//std::cout << "Zadajte bod (X, y) pre zobrazenie ceny opcie (alebo 'exit' pre ukoncenie): " << std::endl;
 //   while (true)
 //   {
	//	std::cin>> std::ws; // ignoruj biele znaky
 //       std::string input;
 //       std::getline(std::cin, input);
 //       if (input == "exit") {
 //           break;
 //       }
 //       double x, y;
 //       try {
 //           size_t pos;
 //           x = std::stod(input, &pos);
 //           y = std::stod(input.substr(pos));
 //       } catch (const std::exception& e) {
 //           std::cerr << "Neplatny vstup. Zadejte dve cisla oddelena mezerou." << std::endl;
 //           continue;
 //       }
 // //      // Najdi indexy i, j pre zadané X a y
 // //      int i = static_cast<int>((x - params.X_l) / solver.getHx());
 // //      int j = static_cast<int>(y / solver.getHy());
 // //      if (i < 0 || i >= Nx || j < 0 || j >= Ny) {
 // //          std::cerr << "Bod mimo rozsah mriežky. Zadejte hodnoty v rozsahu x: [" << params.X_l << ", " << params.X_r << "], y: [0, " << params.y_maxl << "]" << std::endl;
 // //          continue;
 // //      }
	//	//double price = solver.getU()[solver.getMesh().getIndex(i, j)];
	//	//std::cout << "Cena opcie pre x = " << x << ", y = " << y << " je: " << price << std::endl;

	//	// Interpolace ceny pre zadané (x, y)
 //       // Najdi dolný ľavý index
 //       int i = static_cast<int>((x - params.X_l) / solver.getHx());
 //       int j = static_cast<int>(y / solver.getHy());

 //       // Orez na bezpečné hranice (aby i+1 a j+1 existovali)
 //       if (i < 0 || i >= Nx - 1 || j < 0 || j >= Ny - 1) {
 //           std::cerr << "Bod mimo rozsah mriezky..." << std::endl;
 //           continue;
 //       }

 //       // Váhy interpolácie (0.0 až 1.0)
 //       double tx = (x - (params.X_l + i * solver.getHx())) / solver.getHx();
 //       double ty = (y - (j * solver.getHy())) / solver.getHy();

 //       // Štyri rohové hodnoty
 //       const auto& u = solver.getU();
 //       const auto& mesh = solver.getMesh();

 //       double u00 = u[mesh.getIndex(i, j)];
 //       double u10 = u[mesh.getIndex(i + 1, j)];
 //       double u01 = u[mesh.getIndex(i, j + 1)];
 //       double u11 = u[mesh.getIndex(i + 1, j + 1)];

 //       // Bilineárna interpolácia
 //       double price = (1 - tx) * (1 - ty) * u00
 //           + tx * (1 - ty) * u10
 //           + (1 - tx) * ty * u01
 //           + tx * ty * u11;

 //       std::cout << "Cena opcie pre x = " << x << ", y = " << y << " je: " << price << std::endl;
 //   }

      return 0;
}




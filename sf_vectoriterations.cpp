#include "sfnewton.h"



void SFNewton::vector_iterate(double *const x_, int nvars) {
	cout << "Number of iteration variables: " << nvars << endl;
	tick_timer = clock();
	algorithm_ticks = 0;
	gradient_ticks = 0;

	z_method = LIMITED_STORAGE_PSEUDO_HESSIAN; //patrick; 
	//z_method = NEGATIVE_GRADIENT;
	beta_method = CONVERGENCE_DEPENDENT; //patrick's choice
	//beta_method = R_OVER_R_PREV;
	//beta_method = RVDxG_OVER_RVDxPV;
	//beta_method = BETA_ALWAYS1;
	//beta_method = BETA_ALWAYS0;
	//beta_method = Polak_Ribiere; 
	alpha_method = RESET_DEPENDENT; //patrick; 
	//alpha_method = PARABOLA;
	//alpha_method = BISECTIONS; 
	//alpha_method = RR_DIV_PP;

	// x_ is a pointer to an array of iteration values
	// Copy these values into a Vector object
	
	Vector x(1,nvars);
	memcpy(&x[1], x_, sizeof(*x_) * nvars);
	x_addr = x_;

	gradient_call_count = 0;
	int reset_counter = 0;
	
	// x = the iteration variables
	//   = x' + alpha * p
	//   = (the previous values) + (some step length) * (the step direction)
	//
	// p = the step direction
	//   = z + beta * p'
	//   = (a function of the iteration variable x) + (a function of the previous step direction)
	//
	// z = the pseudo Newton step
	//   = -r / A 
	//
	// beta = the conjugate gradient mixing value
	// alpha = the step size factor

	double r_dot_x;
	double r_dot_x_norm;
	//double ys_div_yy; //zit nu in header
	
	// Initialize z to the negative gradient; it's a sensible starting point.
	//		Plus, there's a good chance we do not have enough information yet to call computeZ() properly 
	//		(for instance, we do not have an iteration history yet)
	
	rv = computeGradient(x);
	zv = -rv; 
	pv = zv / (1.0+zv.norm());

	
	alpha = 1.0;
	beta = 0.0;
	
	// Initialize a few more common variables
	Vector x_best = x;
	Vector r_best = rv;
	double r_best_norm = rv.norm();

	int i_prev_imp = 0;
	double r_norm_prev_imp = 0;
	
	cout.precision(4);
	cout.setf(std::ios::scientific);
	
	int i = 0;
	bool continue_iterating = true;
	bool reinitialize_step_direction = false;
	bool reset_to_best_x = false;

	if (rv.norm() < tolerance) {
		cout << "Your starting vector is good enough!" << endl;
		continue_iterating = false;
	}
	
	relative_improvement = 2; // any number > 1.001 suffices for initialization, because we want to start without CG

	while (continue_iterating) {
		bool overflow = false; 
		++i;
		// Sanity check: p must not contain NaN values
		if (isnan(pv.norm())) {
			//cout << endl << "!!! ERROR !!! 	P is not a number!" << endl;
			cout<<"*";
			if (print_common_info) cout << pv << endl;
			//frans
			pv = zv; 
			overflow = true; 
			//break;
		}
		
		// Store the current x and r vectors for later use
		//TODO: Check if this is right (since only the reference to x is copied)
		x_prev = x;
		r_prev = rv;

		alpha = computeAlpha(x);
		if (print_verbose_info) cout << "alpha = " << alpha << endl;
		
		//TODO: Check if this is right (the original x array is lost (which makes the prev TODO work, btw))
		x = x + alpha * pv;
		if (print_verbose_info) cout << "X = " << x << endl;

		rv = computeGradient(x);
		if (print_verbose_info) cout << "R = " << rv << endl;

		// Store the last 'amount_of_stored_vectors' x and r vectors to be used in computeZ()
		if (static_cast<signed int>(rho.size()) >= amount_of_stored_vectors) {
			x_diff.pop_front();
			r_diff.pop_front();
			rho.pop_front();
		}
		x_diff.push_back(x - x_prev);
		r_diff.push_back(rv - r_prev);

		r_dot_x = r_diff.back().dotproduct(r_diff.back());
		cout << "back to back "<< r_dot_x << endl;
		r_dot_x = r_diff.front().dotproduct(r_diff.front());
		cout << "front to front "<< r_dot_x<< endl;
		r_dot_x = x_diff.back().dotproduct(x_diff.back());
		cout << "back to back "<< r_dot_x << endl;
		r_dot_x = x_diff.front().dotproduct(x_diff.front());
		cout << "front to front "<< r_dot_x<< endl;
		
		r_dot_x = r_diff.back().dotproduct(x_diff.back());
		r_dot_x_norm = r_dot_x / r_diff.back().norm() / x_diff.back().norm();
		rho.push_back(1.0 / r_dot_x);
		
		if (print_verbose_info) cout << "|r.x| = " << r_dot_x << endl;
		ys_div_yy= r_dot_x/(r_diff.back().dotproduct(r_diff.back()));
		zv = computeZ(x,ys_div_yy);
		
		if (print_verbose_info) cout << "Z = " << zv << endl;
		
		// Reinitialize every N steps to emulate a full Newton step
		reinitialize_step_direction = (i % nvars == 0); 
		

		// Check various conditions, to see if we need to reinitialize the iteration
		reset_to_best_x = false;
		
		if (isnan(zv.norm())) {
			// The gradient contains NaN values. Reset, otherwise we would have to stop iterating altogether
			reset_to_best_x = true;
		};
		if (r_dot_x_norm < reset_limit_gradient_inversion) {
			// The normalized dotproduct of delta r and delta x is less than a given limit (default -1.0, effectively skips this check)
			reset_to_best_x = true;
			//overflow = true;
			cout << "!"; 
		}
		
		if (reset_to_best_x) {

			// Echo a message according to the problems that arose
			if (print_common_info) {
				cout << "    < Resetting > ";
				if (isnan(zv.norm())) cout << " (Z is not a number) ";
				if (r_dot_x_norm < reset_limit_gradient_inversion) cout << " (|r.x| = " << r_dot_x_norm << " < " << reset_limit_gradient_inversion << ") ";
				cout << endl;
			};
			cout << "<";
			++reset_counter;

			// Clear the storage arrays, because we are changing x to a vector unrelated to recent iterations
			x_diff.clear();
			r_diff.clear();
			rho.clear();
			
			x = x_best;
			rv = r_best;
			zv = computeZ(x,1.0);
			pv = zv;
			
			// Recalculate alpha, and divide it by a factor that depends on the reset_counter, to guarantee a different step after each reset
			alpha = computeAlpha(x) / (reset_counter + 1) / (reset_counter + 1);

			x = x + alpha * pv;
			rv = computeGradient(x);

			

			x_diff.push_back(alpha * pv);
			r_diff.push_back(rv - r_best);
			r_dot_x = r_diff.back().dotproduct(x_diff.back());
			r_dot_x_norm = r_dot_x / r_diff.back().norm() / x_diff.back().norm();
			rho.push_back(1.0 / r_dot_x);
			
			

			zv = computeZ(x,1.0);

			reinitialize_step_direction = true;
		}

		if (reinitialize_step_direction || overflow) {
			overflow = false; 
			pv = zv;
			beta = 0.0;
		}
		else {
			beta = computeBeta(x);
			if (isnan(beta)) {
				if (print_common_info) cout << "    < Resetting > (beta is not a number)" << endl;
				pv = zv;
			}
			else {
				if (print_verbose_info) cout << "beta = " << beta << endl;
				pv = zv + beta * pv;
			}
		}
		if (print_verbose_info) cout << "P = " << pv << endl;
		
		accuracy = rv.norm();
		
		if (accuracy < tolerance) {
			// Stop if we've reached the desired accuracy
			continue_iterating = false;
		}
		if (i >= iterationlimit) {
			// Stop if we've performed the maximum allowed number of iterations
			continue_iterating = false;
		}

		if (rv.norm() < r_best_norm) {

			// Store the best x vector and the associated r vector for two reasons:
			//		- So we have something to fall back on in case of numerical problems (e.g. when z contains NaN values)
			//		- To use as our best estimate if the iteration fails to reach the desired accuracy within the maximum number of allowed iterations

			x_best = x;
			r_best = rv;
			r_best_norm = rv.norm();

			reset_counter = 0;
			if (print_common_info) cout << "    [ Improvement ]\t" << i << "\t" << accuracy << endl;
			
			if (i_prev_imp > 0 && r_best_norm < 1e-2) {
				relative_improvement = pow(r_norm_prev_imp / r_best_norm,  1.0 / (i - i_prev_imp));
			}

			i_prev_imp = i;
			r_norm_prev_imp = r_best_norm;  
		}
		
		if (i % print_iteration_info == 0) {
			cout << " i = " << i << "     best = " << r_best_norm << "     accuracy = " << accuracy << "     alpha = " << alpha << "     beta = " << beta << "     y.s/ys = " << r_dot_x_norm << endl;
		}
		
		if (print_exportable_info || (print_improvement_info && i_prev_imp == i)) {
			cout << i << " " << rv.norm() << " " << alpha << " " << beta << " " << r_dot_x_norm << " ";
			cout << ((reset_to_best_x) ? "reset" : "normal") << " ";
			cout << ((i_prev_imp == i) ? "improvement" : "no");
			cout << endl;
		}
		if (print_verbose_info) cout << endl;
		
	}
	algorithm_ticks += clock() - tick_timer;

	x = x_best;
	rv = r_best;
	
	//if (print_common_info) {
		
		cout << "iteration count = " << i << endl;
		cout << "gradient call count = " << gradient_call_count << endl;
		cout << "ticks spent outside gradient calls = " << algorithm_ticks << endl;
		cout << "ticks spent in gradient calls  = " << gradient_ticks << endl;
		cout << endl << "Newton done." << endl;
	//}
	cout << endl;
}

Vector SFNewton::computeGradient(Vector a) {
	int nvar = a.Length();
	Vector g(1,nvar);
	Vector temp(1,nvar);

	// This method uses memcpy() to be able to calculate gradients of *any* vector, not just x
	// This is needed because the encompassing classes assume the vector to be x
	
	++gradient_call_count;

	memcpy(&temp[1], x_addr, sizeof(*x_addr) * nvar);	// Copy the data held in x into temp
	memcpy(x_addr, &a[1], sizeof(*x_addr) * nvar);		// Copy the data held in a into x

	algorithm_ticks += clock() - tick_timer;
	tick_timer = clock();
	residuals(&g[1], x_addr);									// Put the gradient in the vector part of Vector g
																		// The x_addr is ignored but passed for consistency
	gradient_ticks += clock() - tick_timer;
	tick_timer = clock();

	memcpy(x_addr, &temp[1], sizeof(*x_addr) * nvar);	// Copy the data held in temp back into x
	return g;
}

Vector SFNewton::computeZ(Vector& a, double ys_div_yy) {
	Vector z(1,a.Length());
	deque<double> alphad;
	deque<double> betad;
	
	switch (z_method) {
		
	case NEGATIVE_GRADIENT:
		z = -rv;
		break;
		
	case LIMITED_STORAGE_PSEUDO_HESSIAN:
		
		// Based on "Updating Quasi-Newton Matrices with Limited Storage", Mathematics of Computation, Vol. 35, No. 151 (Jul 1980), pp 773-782, p. 779
		// Note that incr is always 0, due to the way the deques are handled, so the references to mj could technically be replaced
		// by references to mi. It's been kept like this to coincide with the aforementioned article though
		
		
		if (relative_improvement <= 1.001) {
			//cout << "    < warning > relative_improvement <= 1.001..." << endl; 
			z = -rv; 
			break;
		}
		else {
			int incr;
			incr = 0;
			int bound;
			
			if (static_cast<signed int>(rho.size()) <= amount_of_stored_vectors) {
				bound = static_cast<signed int>(rho.size()) - 1;
			}
			else {
				bound = amount_of_stored_vectors - 1;
			}
			if (bound > static_cast<signed int>(rho.size()) - 1) {
				bound = static_cast<signed int>(rho.size()) - 1;
			}
	
			z = rv;
			int mj;
	
			for (int mi = bound; mi >= 0; --mi) {
				mj = mi + incr;
				alphad[mi] = rho[mj] * x_diff[mj].dotproduct(z);
				z = z - alphad[mi] * r_diff[mj];
			}
			z *= ys_div_yy;
			for (int mi = 0; mi <= bound; ++mi) {
				mj = mi + incr;
				betad[mj] = rho[mj] * r_diff[mj].dotproduct(z);
				z = z + x_diff[mj] * (alphad[mi] - betad[mi]);
			}
			z = -z;
			break;
		}
		
	default:
		cout << "    < warning > No Z-method given. Defaulting to -rv..." << endl;
		z = -rv;
		break;
	}
	return z;
};

double SFNewton::computeAlpha(Vector& x) {

	bool verbose_info = false;
	
	switch (alpha_method) {
		
	case ALPHA_ALWAYS1: {
			
			// Always returns 1. Works surprisingly decent, but of course doesn't adapt to ill-conditioned situations
			alpha = 1.0;
			break;
		}
		
	case ALWAYS_INV_ITVARS: {
			
			// Always returns 1 divided by the number of iteration variables N. Works very poorly in general.
			alpha = 1.0 / x.Length();
			break;
		}
	
	case RR_DIV_PP: {
		alpha = rv.dotproduct(rv)/pv.dotproduct(pv);
		if (alpha>1) alpha=1;
		break; 
	}
		
	case PARABOLA: {
			
			// Using the extremum of a parabolic fit to three points performes reasonably well as an alpha method,
			// but has one major drawback: it always uses three gradient calls, which makes this algoritm very time consuming
			// when the number of iteration variables is even somewhat large
			
			// DETERMINING ALPHA BY PARABOLIC FIT
			// 
			// One can fit a parabola to three given points (x, y), using:
			//		x(n) * x(n) * a + x(n) * b + c = y(n)
			//
			// For x1 = 0, x2 = 1/2 and x3 = 1, these equations reduce to:
			//		              c = y1
			//		1/4a + 1/2b + c = y2
			//		   a +    b + c = y3
			//
			// Solving for a and b gives:
			//		b = 4y2 - y3 - 3y1
			//		a = y3 - y1 - b
			//
			// The extremum x(extr) of a parabola can be found at:
			//		x(extr) = -b / 2a
			//
			// Substitution yields:
			//		x(extr) = (3y1 - 4y2 + y3) / (4y1 - 8y2 + 4y3)
			// 
			// If a is equal to or less than zero, the parabola has a maximum instead of a minimum.
			// In that case, 1 / nvars is returned for alpha
			
			double y1 = computeGradient(x).norm();
			double y2 = computeGradient(x + pv / 2).norm();
			double y3 = computeGradient(x + pv).norm();
			double b = 4 * y2 - y3 - 3 * y1;
			double a = y3 - y1 - b;
			
			double x_peak = -b / a / 2;
			
			if (a < 0 && x_peak < 0) alpha = 1;
			if (a < 0 && x_peak >= 0 && x_peak <= 1) alpha = (y1 > y3) ? 1.0 : 1.0 / x.Length();
			if (a < 0 && x_peak > 1) alpha = 1.0 / x.Length();
			
			if (a > 0 && x_peak < 0) alpha = 1.0 / x.Length();
			if (a > 0 && x_peak >= 0 && x_peak <= 1) alpha = x_peak;
			if (a > 0 && x_peak > 1) alpha = 1;
			
			if (verbose_info) cout << "(y1, y2, y3)  = " << y1 << ", " << y2 << ", " << y3 << endl;
			if (verbose_info) cout << "(a, b, alpha) = " << a << ", " << b << ", " << alpha << endl;
			break;
		}
		
	case R_DEPENDENT: {
			
			// This method attempts to adapt alpha dynamically to the current circumstances:
			//		- Increase alpha when things seem to be going smoothly
			//		- Decrease alpha when trouble arises
			//	However, it takes a lot of parameters, which are not very well optimized through experimentation
			// As it stands now, it reduces alpha ruthlessly, and very small alpha values (< 1e-6) lead to very small steps
			// This method can potentially bring the main algoritm to a grinding halt.
			
			double alpha_residue_crit = 2;
			double alpha_dec_factor = 2;
			double alpha_reward_factor = 1.1;
			
			double r_prev_norm, r_new_norm;

			// Check if the |r(new)| is not substantially larger than |r(old)|
			r_prev_norm = rv.norm();
			r_new_norm = computeGradient(x + pv * alpha).norm();

			if (alpha < numeric_limits<double>::epsilon() || isnan(alpha)) {
				alpha = 1.0;
			}
			if (r_new_norm > r_prev_norm * alpha_residue_crit) {
				alpha /= alpha_dec_factor;
				if (verbose_info) cout << "    - decreased alpha to " << alpha << " because " << r_new_norm << " >> " << r_prev_norm << " -" << endl;
			}
			else if (r_new_norm < r_prev_norm) {
				alpha *= alpha_reward_factor;
				if (alpha > 1) alpha = 1;
				if (verbose_info) cout << "    - increased alpha to " << alpha << " -" << endl;
			}
			break;
		}
		
	case BISECTIONS: {
			double alpha_l = 0.0;
			double alpha_r = 1.0;

			double new_alpha = 0.0;
			double norm_l = 0.0;
			double norm_r = 0.0;
			
			norm_l = computeGradient(x + pv * alpha_l).norm();
			norm_r = computeGradient(x + pv * alpha_r).norm();
			for (int bicount = 0; bicount < 4; ++bicount) {
				if (verbose_info) cout << "(L,R) (" << alpha_l << ", " << alpha_r << ") = " << norm_l << ", " << norm_r << endl;
				new_alpha = (alpha_l + alpha_r) / 2;

				if (norm_l > norm_r) {
					alpha_l = new_alpha;
					norm_l = computeGradient(x + pv * alpha_l).norm();
				}
				else {
					alpha_r = new_alpha;
					norm_r = computeGradient(x + pv * alpha_r).norm();
				}
			}
			if (verbose_info) cout << "(L,R) (" << alpha_l << ", " << alpha_r << ") = " << norm_l << ", " << norm_r << endl;
			
			// Choose a final alpha to return
			if (norm_l > norm_r) {
				alpha = alpha_r;
			}
			else {
				alpha = alpha_l;
			}
			if (alpha < numeric_limits<double>::epsilon()) {
				alpha = (alpha_l + alpha_r) / 20;
			}
			break;
		}
		
	case RESET_DEPENDENT: {
			
			// This method does not compute a whole new alpha value, but just passes the current one back, after bringing it closer to 1.0
			// Works well with the alpha reduction that occurs on resets (hence the name), although the current divisor 5 might be somewhat low:
			// a higher value would mean more cautious steps after a reset
			
			alpha += (1.0 - alpha) / 5;
			break;
		}

	default: {
			cout << "    < warning > No alpha-method given. Defaulting to 1.0 ..." << endl;
			alpha = 1.0;
			break;
		}
	}

	// Sanity checks on alpha
	if (alpha < numeric_limits<double>::epsilon() || isnan(alpha)) {
		// If alpha is not a number, less than or (almost) equal to zero, return 1 divided by the number of iteration variables, an arbitrary small number
		alpha = 1.0 / x.Length();
	}

	if (alpha > 1.0) {
		// The system is not well defined for alpha > 1
		alpha = 1.0;
	}
	
	return alpha;
}

double SFNewton::computeBeta(Vector& x) {
	switch (beta_method) {
		
	case BETA_ALWAYS1: {
	
			// There is no theoretical basis for this beta method ...
			beta = 1.0;
			break;
		}
	case BETA_ALWAYS0: {
	
			// A beta of zero basically strips the algoritm of its conjugate gradient component, and turns it into a (modified) steepest descent method
			beta = 0.0;
			break;
		}
		
	case RVDxG_OVER_RVDxPV: {
	
			// Based on "Conjugate Gradient Methods with Inexact Searches", D. Shanno, Mathematics of Operations Research, Vol. 3, No. 3, Aug. 1978, 
			//     pp 244-256, p244
			
			double rd_dot_r = (rv - r_prev).dotproduct(rv);
			double rd_dot_pv = (rv - r_prev).dotproduct(pv);
			beta = rd_dot_r / rd_dot_pv;
			break;
		}
 
	case Polak_Ribiere: {
			double rd_dot_r = (rv - r_prev).dotproduct(rv);
			double r_dot_r = r_prev.dotproduct(r_prev); 
			beta = rd_dot_r/r_dot_r;
			if (beta < 0) beta = 0 ;  
			//if (beta > 100) beta = 0;
			break;
		}
		
	case R_OVER_R_PREV: {
	
			// Based on "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", J.R. Shewchuk, Edition 1 1/4, August 4, 1994, 
			//     http://www.cs.cmu.edu/~jrs/jrspapers.html
			
			beta = rv.norm() / r_prev.norm();
			
			break;
		}
		
	case CONVERGENCE_DEPENDENT: {
		
			// A combination of BETA_ALWAYS0 (no conjugate gradients) and RVDxG_OVER_RVDxPV (conjugate gradients)
			// It appears that CG works well when we are *very* close to the minimum, but not when we are either meandering near the starting point
			// or when the routine is in a steep path towards the minimum
			// The relative_improvement is a measure of the strength of the most recent improvement. If this value is less than 1.001, it's very likely
			// that we are in a plateau near the minimum, so we should use CG
			
			if (relative_improvement > 1.001) {
			//if (accuracy > tolerance*100) {
				beta = 0;
			}
			else {
				double rd_dot_r = (rv - r_prev).dotproduct(rv);
				double rd_dot_pv = (rv - r_prev).dotproduct(pv);
				beta = rd_dot_r / rd_dot_pv;
				//beta = rd_dot_pv / rd_dot_r;
			}
			break;
		}	

	default: {
			cout << "    < warning > No beta-method given. Defaulting to 0.0 ..." << endl;
			beta = 0.0;
			break;
		}
	};
	
	return beta;
}


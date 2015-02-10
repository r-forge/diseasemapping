cholmod_common c;

SEXP gmrfLik(SEXPR QR, SEXP QcholR, SEXP obsCovR,  SEXP propNuggetR)
{
    CHM_FR Qchol = AS_CHM_FR(QcholR), Vchol;
    CHM_SP Q = AS_CHM_SP__(QR);


    CHM_SP B = AS_CHM_SP__(obsCovR);
    int sys = asInteger(system);

		double *propNugget = REAL(propNuggetR);
		int Nnugget = length(propNuggetR);

    R_CheckStack();

    Vchol = cholmod_copy_factor(Qchol, &c);

    Vchol = chm_factor_update(Vchol, Q, propNugget[0]);

	// Vx = solve(Vchol, X,system='L')
	Vx = cholmod_solve();


		
}

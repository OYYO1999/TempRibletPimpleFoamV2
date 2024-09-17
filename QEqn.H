// Compute Q for Esource
// now

scalar Tao_relax = 0.3 * delta.value() / u_tao.value();
scalar Tmean = 0.0;
scalar rhomean = 0.0;
scalar Cvmean = 0.0;

scalar V = gSum(mesh.V());

scalarField cellVolume = mesh.V();

volScalarField T_q = thermo.T();
volScalarField rho_q = thermo.rho();
volScalarField Cv_q = thermo.Cv();

for(label celli=0; celli < mesh.C().size(); celli++)
{
	scalar volCell = cellVolume[celli];
	Tmean += T_q[celli] * volCell;
	rhomean += rho_q[celli] * volCell;
	Cvmean += Cv_q[celli] * volCell;
}

Tmean /= V;
rhomean /= V;
Cvmean /= V;

reduce(Tmean, sumOp<scalar>());//需要全部规约才行，并行才不会出问题
reduce(rhomean, sumOp<scalar>());//
reduce(Cvmean, sumOp<scalar>());//
reduce(V, sumOp<scalar>());//

Q.primitiveFieldRef() = rhomean * Cvmean * (TRef.value() - Tmean) / Tao_relax;

define MAX_ITR;
define kappa_param 
define tension_param 
MAX_ITR := 20;
kappa_param := 25
tension_param := 0.25

define h_zero_min;
define h_zero_max;
define lambda_min; 
define lambda_max;
h_zero_min := -0.3;
h_zero_max := 0.6;
lambda_min := 1; 
lambda_max := 50;


read "obj.cmd"

mean_curv.modulus := kappa_param;
default_area.modulus := tension_param;
recalc;
print "Total Energy:";
print total_energy;
show_all_quantities;
Q;
QUIETGO;

procedure edge_fixing()
{
    // Step 1: Fix edges connected to fixed vertices
    foreach vertex vv where vv.fixed do
    {
        foreach edge ee where ee.vertex[1].id == vv.id or ee.vertex[2].id == vv.id do
        {
            fix ee;  // Fix the edge
        };
    };
    // Step 2: Fix vertices connected to fixed edges
    foreach edge ee where ee.fixed do
    {
        fix ee.vertex[1];  // Fix the first vertex of the edge
        fix ee.vertex[2];  // Fix the second vertex of the edge
    };
    set edge color blue where fixed;
}


procedure unfixing()
{
	foreach vertex vv where vv.fixed do
	{
		unfix vv
	};
	foreach edge ee where ee.fixed do
	{
		unfix ee
	};
}


procedure toggle_fixation()
{
    // Step 1: Toggle fixation for all fixed edges
    foreach edge ee do
    {
        if (ee.fixed) then
            unfix ee // Unfix the edge
        else
            fix ee;    // Fix the edge
    };

    // Step 2: Toggle fixation for all fixed vertices
    foreach vertex vv do
    {
        if (vv.fixed) then
            unfix vv  // Unfix the vertex
        else
            fix vv;    // Fix the vertex
    };
	set edge color black;
    set edge color blue where fixed;
}

procedure smooth_rim()
{
	toggle_fixation();
	default_area.modulus := 0;
	g100;
	toggle_fixation();
}

procedure set_sigma(integer lambda)
{
	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;
}

gogo_hessian :=
{
	g 40;
	V 20; u 4; V 20; u 2;
	g 10;
	hessian_seek; 
	hessian_seek;
	g 30;
}

gogo_h := {
	gogo_hessian;
	gogo_hessian;
	gogo_hessian;
	gogo_hessian;
	gogo_hessian;
	gogo_hessian;
}

procedure mini_hessian()
{   
	local E_EPS;
	local old_energy;
	E_EPS := 1e-2;
	conj_grad ON;
	old_energy := total_energy;
	for( it := 1; it<MAX_ITR; it ++)
	{	
		printf "Iteration Number %d \n", it;
		t 0.5; u 2;
		gogo_hessian;
        if old_energy < total_energy then
		{
			gogo_hessian;
            printf "BREAK: old_energy: %.3f  < total_energy: %.3f \n", old_energy, total_energy;
			break;
		};
		if abs(total_energy-old_energy)< E_EPS then
		{
            print "BREAK: abs(total_energy-old_energy) < E_EPS \n";
			break;
		};
		if scale < 1e-12 then
		{
            print "BREAK: scale < 1e-12 \n";
			break;
		};
		old_energy := total_energy;
        l 1;
		K 1e-2;
		t 1e-3; u 2;
	};
    l 0.5;
    K 1e-2;
	t 1e-3; u 2;
    V 5; u 2;
	g20;
	hessian_seek;
	hessian_seek;
	g50;
}

procedure create_obj(integer h_it, integer lam_it, integer grid_sice)
{	
	print "start";
	edge_fixing();
	print "Edges fixed";
	smooth_rim();
	print "Rim smoothed";
	h_zero := h_zero_min + h_it * (h_zero_max - h_zero_min) / (grid_sice - 1);
	
	ex_l := lam_it / (grid_sice - 1); 
	lambda := lambda_min * pow(lambda_max / lambda_min, ex_l);


	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;

	mini_hessian();

	
	printf "Total Energy= %.3f \n h_zero= %.3f \n kappa= %.3f \n sigma= %.3f \n",total_energy, h_zero ,mean_curv.modulus, default_area.modulus;
	Q;
	obj >>> sprintf "Min_Necks_temp/Neck%d_%d%d.obj",Neck_nr,h_it,lam_it;
		

	
}

procedure create_obj_test(integer h_it, integer lam_it, integer grid_sice)
{	
	print "start";
	print "Edges fixed";
	h_zero := h_zero_min + h_it * (h_zero_max - h_zero_min) / (grid_sice - 1);
	
	ex_l := lam_it / (grid_sice - 1); 
	lambda := lambda_min * pow(lambda_max / lambda_min, ex_l);


	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;

	mini_hessian();

	
	printf "Total Energy= %.3f \n h_zero= %.3f \n kappa= %.3f \n sigma= %.3f \n",total_energy, h_zero ,mean_curv.modulus, default_area.modulus;
	Q;
	obj >>> sprintf "Min_Necks_temp/Neck%d_%d%d.obj",Neck_nr,h_it,lam_it;
		

	
}


procedure create_obj_lambda(integer h_it, integer lambda, integer grid_sice)
{	
	print "start";
	print "Edges fixed";
	print "Rim smoothed";
	h_zero := h_zero_min + h_it * (h_zero_max - h_zero_min) / (grid_sice - 1);



	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;

	mini_hessian();

	
	printf "Total Energy= %.3f \n h_zero= %.3f \n kappa= %.3f \n sigma= %.3f \n",total_energy, h_zero ,mean_curv.modulus, default_area.modulus;
	Q;
	obj >>> sprintf "Min_Necks/Min_Necks_temp/Neck%d_%d_%d.obj",Neck_nr,lambda,h_it;
		

	
}

procedure create_obj_const_lam(integer h_it, integer lambda, integer grid_sice, integer interval_number)
{	
	if interval_number == 0 then 
	{
		h_zero_min := -0.1;
		h_zero_max := 0.4;
	};
	if interval_number == 1 then 
	{
		h_zero_min := -0.2;
		h_zero_max := 0.2;
	};
	if interval_number == 2 then 
	{
		h_zero_min := -0.1;
		h_zero_max := 0.1;
	};
	if interval_number == 3 then 
	{
		h_zero_min := -0.025;
		h_zero_max := 0.012;
	};
	if interval_number == 4 then 
	{
		h_zero_min := -0.075;
		h_zero_max := 0.012;
	};


	h_zero := h_zero_min + h_it * (h_zero_max - h_zero_min) / (grid_sice - 1);
	

	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;
	

	mini_hessian();

	
	printf "Total Energy= %.3f \n h_zero= %.3f \n kappa= %.3f \n sigma= %.3f \n",total_energy, h_zero ,mean_curv.modulus, default_area.modulus;
	Q;

	obj >>> sprintf "Min_Necks/Min_Necks_temp/Neck%d_%d_%d_%d.obj",Neck_nr,lambda,h_it,interval_number;
	
}





procedure iterate_h_zero_hessian_seek_obj(integer grid_sice, integer lambda)
{	
	sigma := kappa_param/(lambda**2); 
	default_area.modulus := sigma;
	
	for (h_it := 0; h_it < grid_sice; h_it ++)
	{	
		REPLACE_LOAD sprintf "code_Necks_fe/code_Neck%d.fe", Neck_nr;
		edge_fixing();

		h_zero := h_zero_min + h_it * (h_zero_max - h_zero_min) / (grid_sice - 1);

		mini_hessian();
		printf "Total Energy= %.3f \n h_zero= %.3f \n kappa= %.3f \n sigma= %.3f \n",total_energy, h_zero ,mean_curv.modulus, default_area.modulus;
		Q;
		obj >>> sprintf "Min_Necks/Neck%d/Neck%d_%d_%d.obj",Neck_nr,Neck_nr,lambda,h_it;
		

	};

}


##Rstan code for pre and post 1980 model.extract()

stanmodelcode_0813_Q_PIC_power_prepost = '

data {  
  int nl;                     		//3600 subwatersheds separated by year       
  int nr;   				//483 incremental watershed-years
  int nd;                			//985 dischargers associated with "nr"
  vector [nr] load;			//WRTDS load at LMS
  vector [nr] SD;			//SD calculated for WRTDS load
  vector [nl] chick;			//chickens in subwatershed
  vector [nl] cow;			//cow in subwatershed
  vector [nl] hog;			//hogs (swine) in subwatershed
  int wsd [nr];			//count variable for LMSs
  int wshed_size;			//# of LMS watersheds
    vector [nr] increm_area;  	//Incremental area for each loading station
  vector [nr] av_prec;		//normalized precipitation
  vector [nr] av_prec2;		//scaled precipitation
  vector [nr] up_t_load1;		//upstream TN loading for nested wsds
  vector [nr] up_t_load2;		//upstream TN loading for nested wsds
  vector [nr] up_t_load3;		//upstream TN loading for nested wsds
  vector [nr] str_loss_load1;	//stream losses of upstream TN loading
  vector [nr] str_loss_load2;	//stream losses of upstream TN loading
  vector [nr] str_loss_load3;	//stream losses of upstream TN loading
  vector [nr] res_loss_load1;	//reservoir losses of upstream TN loading
  vector [nr] res_loss_load2;	//reservoir losses of upstream TN loading
  vector [nr] res_loss_load3;	//reservoir losses of upstream TN loading
  vector [nd] d_loss_str;		//stream losses of dischargers to LMSs
  vector [nd] d_loss_res;		//reservoir losses of dischargers to LMSs
  vector [nd] d_vals;		//dischargers TN loadings
  vector [nl] l_loss_str;		//stream losses for subwatersheds
  vector [nl] l_loss_res;		//reservoir losses for subwatersheds
  vector [nl] ag;			//agriculture in subwatershed
  vector [nl] devpre;		//pre-1980 urban in subwatershed
  vector [nl] devpost;		//post-1980 urban in subwatershed
  vector [nl] wild;			//undeveloped in subwatershed
  vector [nl] tot_l;			//total size of subwatershed
  int l_start [nr];			//count variables to link subwatersheds to LMSs
  int l_end [nr];			//count variables to link subwatersheds to LMSs
  int d_start [nr];			//count variables to link dischargers to LMSs
  int d_end [nr];			//count variables to link dischargers to LMSs
}

transformed data{
}

parameters {  
  real<lower =0> Be_a;        			//Agriculture export
  real<lower =0> Be_d_pre;        			//pre-1980 developed export
  real<lower =0> Be_d_post;        			//post-1980 developed export
  real<lower =0> Be_w;       			//Undeveloped export
  real<lower =0, upper = 1> Be_ch;       		//Chickens export
  real<lower =0, upper = 1> Be_h;      		//Hogs (swine)
  real<lower =0, upper = 1> Be_cw;       		//Cows
  real <lower =0, upper = 1> Sn;			//Stream retention
  real <lower =1, upper = 60> Sn2;		//Reservoir retention
  real <lower =0, upper = 0.40> PIC_q;		//PIC for stream flow
  vector<lower =0, upper = 10> [7] pic_p;		//PIC for land classes
  real <lower=0> Be_dch;      			//Point source discharge coefficient 
  real<lower=0, upper = 2> sigma_res;		//Model residual sigma
  real <lower=0, upper = 5> sigma_w;		//Random effect sigma
  vector [wshed_size] alpha;			//# of watersheds
    real<lower = 0, upper = 2> sigma_B1;		//PIC sigma
  real<lower = 0, upper = 3> Bp_mean;		//PIC mean
  vector [nr] ly;					//log unknown true loads
}

transformed parameters {
}

model {
  vector [nr] tot; 			// Sum of all loadings from all sources
  vector [nr] sigma;		//sigma for watershed random effects
  vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
  vector [nr] A;			//To compile Agriculture with PIC
  vector [nr] Dpr;			//To compile pre-1980 urban with PIC
  vector [nr] Dpt;			//To compile post-1980 urban with PIC
  vector [nr] W;			//To compile undeveloped urban with PIC
  vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
  vector [nr] tot_loss;
  int w;
  vector [nr] alpha_vals;		// Watershed indicator
  real t;
  
  vector [nr] A;			//Agriculture vector
  vector [nr] D_lc_pre;		//pre-1980  vector		
  vector [nr] D_lc_post;		//post-1980  vector
  vector [nr] W_lc;		//Undeveloped  vector
  vector [nr] Disch;		//point source dischargers
  vector [nr] C;   			//chickens for adding PIC
  vector [nr] H;			//hogs for adding PIC
  vector [nr] Cw;			//cows for adding PIC
  vector [nr] C_r;   		//chickens for aggregating subwatersheds
  vector [nr] H_r;			//swine for aggregating subwatersheds
  vector [nr] Cw_r;		//cows  for aggregating subwatersheds
  
  vector [nr] ly_hat;			// log of TN load with transformation and offset
  vector [nr] y;				//TN load
  
  // Loop to determine export for each watershed-year
  for(i in 1:nr){
    
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
    A_lc[i]= sum((ag[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    D_lc_pre[i]= sum((devpre[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    D_lc_post[i]= sum((devpost[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
    C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
    
    //Adding precipitation impact coefficient to land export 
    A[i]=  Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc[i];
    Dpr[i] = Be_d_pre * pow(av_prec2[i],pic_p[2]) .* D_lc_pre[i]; 
    Dpt[i] = Be_d_post * pow(av_prec2[i],pic_p[3]) .* D_lc_post[i];
    D[i]   = Dpr[i]+Dpt[i];
    W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];
    C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
    H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
    Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];
    Dch[i] = Be_dch * Disch[i];
  }
  
  //Loop to determine random effect for each watershed
  for (i in 1:nr){
    w= wsd[i];
    sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
    alpha_vals[i] = alpha[w];
  }
  
  //Sum loadings from all sources
  tot =  A + D + W + C + H + Cw + Dch; 
  
  //Add random effects to source loadings and subtract losses from upstream loads
  y_hat = tot + 10000 * alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));
  
  //Priors for the model
  Be_a ~ normal(900,700);  			//Prior for agriculture, TN/km2/yr
  Be_d_pre ~ normal(800,300);  			//Prior for pre-1980 development, TN/km2/yr
  Be_d_post ~ normal(800,300);  			//Prior for post-1980 development, TN/km2/yr
  Be_w ~ normal(200,200);   			//Prior for undeveloped, TN/km2/yr
  Be_ch ~ normal(0.001,0.0003);  			//Prior for chickens (TN/count)
  Be_h ~ normal(0.04,0.02); 		 	//Prior for swine  (TN/count)
  Be_cw ~ normal(0,5);  				//Uninformed Prior for cow
  Be_dch ~ normal(1,.10);   			//Prior for point source delivery, unitless
  sigma_res ~ normal(0,20);   			//st error of the model
  sigma_w ~ normal(0,300);     			//st. deviation of random effect hyperdistribution
  alpha ~ normal(0,sigma_w);			//watershed random effects
  sigma_B1 ~ normal(0,1);  			//st. deviation for PIC hyperdistribution
  Bp_mean ~ normal(1,1);  			//mean for PIC hyperdistribution
  pic_p[1] ~ normal(Bp_mean,sigma_B1);  	//PIC ag distribution
  pic_p[2] ~ normal(Bp_mean,sigma_B1);  	//PIC pre distribution
  pic_p[3] ~ normal(Bp_mean,sigma_B1);  	//PIC post distribution
  pic_p[4] ~ normal(Bp_mean,sigma_B1);  	//PIC undeveloped distribution
  pic_p[5] ~ normal(Bp_mean,sigma_B1);  	//PIC chicken distribution
  pic_p[6] ~ normal(Bp_mean,sigma_B1);  	//PIC hog distribution
  pic_p[7] ~ normal(Bp_mean,sigma_B1);  	//PIC cow distribution
  Sn ~ normal(.14,.05);				//Prior for stream retention rate, 1/d
  Sn2 ~ normal(11,2);				//Prior for waterbody retention rate, m/y
  PIC_q ~ normal(0,1);				//PIC for retention
  
  // with log transformation 
  ly_hat=log((y_hat ./10000) + 10);   	// vector to get log of TN load (y_hat)
  ly ~ normal(ly_hat,sigma_res);        	//parameter that calibrates ly_hat (log of TN load) with ly () 
  y=exp(ly)-10;               			// parameter that calibrates to load
  load ~ normal(y,SD);               	 	// load = WRTDS estimate, SD = WRTDS sd
}

generated quantities {

}
'




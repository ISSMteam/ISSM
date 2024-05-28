function paterson(temperature){
//PATERSON - figure out the rigidity of ice for a given temperature
//
//   rigidity (in s^(1/3)Pa) is the flow law paramter in the flow law sigma=B*e(1/3) (Paterson, p97). 
//   temperature is in Kelvin degrees
//
//   Usage:
//      var rigidity=paterson(temperature)

	//variables:
	var T=[];

	if (ArrayAnyBelowStrict(temperature,0)){
		throw Error('input temperature should be in Kelvin (positive)');
	}
	
	T=temperature.slice(0);
	for(var i=0;i<temperature.length;i++)T[i]=temperature[i]-273.15;

	//The routine below is equivalent to:

	// n=3; T=temperature-273;
	// //From paterson,
	// Temp=[0;-2;-5;-10;-15;-20;-25;-30;-35;-40;-45;-50];
	// A=[6.8*10^-15;2.4*10^-15;1.6*10^-15;4.9*10^-16;2.9*10^-16;1.7*10^-16;9.4*
	// 10^-17;5.1*10^-17;2.7*10^-17;1.4*10^-17;7.3*10^-18;3.6*10^-18];;//s-1(kPa-3)
	// //Convert into rigidity B
	// B=A.^(-1/n)*10^3; //s^(1/3)Pa
	// //Now, do a cubic fit between Temp and B: 
	// fittedmodel=fit(Temp,B,'cubicspline');
	// rigidity=fittedmodel(temperature);

	var rigidity=NewArrayFill(T.length,0);
	
	for (var i=0;i<T.length;i++){
		
		if(T[i]<=-45)              rigidity[i]=Math.pow(10,8)*(-0.000292866376675*Math.pow(T[i]+50,3)+ 0.011672640664130*Math.pow(T[i]+50,2)  -0.325004442485481*(T[i]+50)+  6.524779401948101);
		if(-45<=T[i] & T[i]<-40)   rigidity[i]=Math.pow(10,8)*(-0.000292866376675*Math.pow(T[i]+45,3)+ 0.007279645014004*Math.pow(T[i]+45,2)  -0.230243014094813*(T[i]+45)+  5.154964909039554);
		if(-40<=T[i] & T[i]<-35)   rigidity[i]=Math.pow(10,8)*(0.000072737147457*Math.pow(T[i]+40,3)+  0.002886649363879*Math.pow(T[i]+40,2)  -0.179411542205399*(T[i]+40)+  4.149132666831214);
		if(-35<=T[i] & T[i]<-30)   rigidity[i]=Math.pow(10,8)*(-0.000086144770023*Math.pow(T[i]+35,3)+ 0.003977706575736*Math.pow(T[i]+35,2)  -0.145089762507325*(T[i]+35)+  3.333333333333331);
		if(-30<=T[i] & T[i]<-25)   rigidity[i]=Math.pow(10,8)*(-0.000043984685769*Math.pow(T[i]+30,3)+ 0.002685535025386*Math.pow(T[i]+30,2)  -0.111773554501713*(T[i]+30)+  2.696559088937191);
		if(-25<=T[i] & T[i]<-20)   rigidity[i]=Math.pow(10,8)*(-0.000029799523463*Math.pow(T[i]+25,3)+ 0.002025764738854*Math.pow(T[i]+25,2)  -0.088217055680511*(T[i]+25)+  2.199331606342181);
		if(-20<=T[i] & T[i]<-15)   rigidity[i]=Math.pow(10,8)*(0.000136920904777*Math.pow(T[i]+20,3)+  0.001578771886910*Math.pow(T[i]+20,2)  -0.070194372551690*(T[i]+20)+  1.805165505978111);
		if(-15<=T[i] & T[i]<-10)   rigidity[i]=Math.pow(10,8)*(-0.000899763781026*Math.pow(T[i]+15,3)+ 0.003632585458564*Math.pow(T[i]+15,2)  -0.044137585824322*(T[i]+15)+  1.510778053489523);
		if(-10<=T[i] & T[i]<-5)    rigidity[i]=Math.pow(10,8)*(0.001676964325070*Math.pow(T[i]+10,3)-  0.009863871256831*Math.pow(T[i]+10,2)  -0.075294014815659*(T[i]+10)+  1.268434288203714);
		if(-5<=T[i] & T[i]<-2)     rigidity[i]=Math.pow(10,8)*(-0.003748937622487*Math.pow(T[i]+5,3)+0.015290593619213*Math.pow(T[i]+5,2)  -0.048160403003748*(T[i]+5)+  0.854987973338348);
		if(-2<=T[i])              rigidity[i]=Math.pow(10,8)*(-0.003748937622488*Math.pow(T[i]+2,3)-0.018449844983174*Math.pow(T[i]+2,2)  -0.057638157095631*(T[i]+2)+  0.746900791092860);

		//Now make sure that rigidity is positive
		if(rigidity[i]<0)          rigidity[i]=Math.pow(10,6);
	}
	return rigidity;
}

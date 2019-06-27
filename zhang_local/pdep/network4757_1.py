species(
    label = 'C#CC([O])O[C]=C(22460)',
    structure = SMILES('C#CC([O])O[C]=C'),
    E0 = (348.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,368.182,368.184,368.187,368.191],'cm^-1')),
        HinderedRotor(inertia=(0.177281,'amu*angstrom^2'), symmetry=1, barrier=(17.0532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752112,'amu*angstrom^2'), symmetry=1, barrier=(72.3524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251165,'amu*angstrom^2'), symmetry=1, barrier=(24.1614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09411,0.0642325,-7.75358e-05,4.88795e-08,-1.21868e-11,42022.7,25.9343], Tmin=(100,'K'), Tmax=(982.249,'K')), NASAPolynomial(coeffs=[12.3523,0.0183861,-7.52344e-06,1.36112e-09,-9.2552e-14,39811.1,-28.1781], Tmin=(982.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(CCOJ)"""),
)

species(
    label = 'C=C=O(598)',
    structure = SMILES('C=C=O'),
    E0 = (-59.3981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52746,0.00708371,9.17709e-06,-1.64254e-08,6.71115e-12,-7123.94,5.7438], Tmin=(100,'K'), Tmax=(956.683,'K')), NASAPolynomial(coeffs=[5.76495,0.00596559,-1.98486e-06,3.52744e-10,-2.51619e-14,-7929,-6.92178], Tmin=(956.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.3981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-(Cdd-O2d)HH)"""),
)

species(
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C1OC1O[C]=C(23834)',
    structure = SMILES('[CH]=C1OC1O[C]=C'),
    E0 = (375.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801595,0.0520979,-7.91119e-06,-4.81493e-08,2.80591e-11,45266.3,22.7651], Tmin=(100,'K'), Tmax=(922.121,'K')), NASAPolynomial(coeffs=[24.4337,-0.00282406,4.01532e-06,-8.03633e-10,4.92164e-14,38884.6,-100.301], Tmin=(922.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C(=C)OC1[O](23835)',
    structure = SMILES('[CH]=C1C(=C)OC1[O]'),
    E0 = (308.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.308,0.0512822,-3.84497e-05,1.13154e-08,-4.34523e-13,37250.1,19.6504], Tmin=(100,'K'), Tmax=(1136.92,'K')), NASAPolynomial(coeffs=[14.2435,0.0166528,-7.11752e-06,1.36118e-09,-9.67813e-14,33605.5,-47.5086], Tmin=(1136.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC(=O)O[C]=C(23836)',
    structure = SMILES('C#CC(=O)O[C]=C'),
    E0 = (245.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,2175,525,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21009,0.0685364,-0.000116718,1.02769e-07,-3.45134e-11,29567.1,24.9358], Tmin=(100,'K'), Tmax=(879.227,'K')), NASAPolynomial(coeffs=[8.15241,0.0214739,-1.002e-05,1.84289e-09,-1.22561e-13,28944.6,-4.26094], Tmin=(879.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'C#COC([O])C#C(23837)',
    structure = SMILES('C#COC([O])C#C'),
    E0 = (318.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,756.667,763.333,770,3350,3450,2000,2200,1380,1390,370,380,2900,435,2100,2250,500,550,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58055,'amu*angstrom^2'), symmetry=1, barrier=(36.34,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58108,'amu*angstrom^2'), symmetry=1, barrier=(36.3521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58119,'amu*angstrom^2'), symmetry=1, barrier=(36.3546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844716,0.0686034,-8.85762e-05,5.63018e-08,-1.38447e-11,38374.3,21.1558], Tmin=(100,'K'), Tmax=(1004.61,'K')), NASAPolynomial(coeffs=[14.8811,0.0127138,-5.12415e-06,9.20816e-10,-6.25553e-14,35554.1,-46.6258], Tmin=(1004.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Ct-CtCs) + group(Ct-CtOs) + group(Ct-CtH) + group(Ct-CtH) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = 'C=[C]OC=O(4555)',
    structure = SMILES('C=[C]OC=O'),
    E0 = (-54.4552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,180,877.069],'cm^-1')),
        HinderedRotor(inertia=(1.02066,'amu*angstrom^2'), symmetry=1, barrier=(23.467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344328,'amu*angstrom^2'), symmetry=1, barrier=(23.4666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43675,0.0295367,-1.86127e-05,3.82274e-09,2.5975e-13,-6489.17,17.8736], Tmin=(100,'K'), Tmax=(1220.16,'K')), NASAPolynomial(coeffs=[9.68768,0.0123156,-5.4936e-06,1.05394e-09,-7.43104e-14,-8746.16,-20.5486], Tmin=(1220.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.4552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CJO)"""),
)

species(
    label = 'C#C[C](O)O[C]=C(23838)',
    structure = SMILES('C#C[C](O)O[C]=C'),
    E0 = (328.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687593,0.0687235,-7.94751e-05,4.13407e-08,-6.96812e-12,39580.8,27.2175], Tmin=(100,'K'), Tmax=(902.127,'K')), NASAPolynomial(coeffs=[16.9911,0.00940524,-2.4114e-06,3.2905e-10,-1.96132e-14,36111.5,-52.6835], Tmin=(902.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([O])C#C(23839)',
    structure = SMILES('[CH]=COC([O])C#C'),
    E0 = (355.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32657,'amu*angstrom^2'), symmetry=1, barrier=(30.5004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32761,'amu*angstrom^2'), symmetry=1, barrier=(30.5243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3288,'amu*angstrom^2'), symmetry=1, barrier=(30.5517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672563,0.0676852,-7.34966e-05,3.58969e-08,-5.64696e-12,42927.4,24.7175], Tmin=(100,'K'), Tmax=(945.383,'K')), NASAPolynomial(coeffs=[17.1149,0.0107662,-3.25651e-06,5.18653e-10,-3.42819e-14,39253.2,-56.6736], Tmin=(945.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=C([O])OC=C(23840)',
    structure = SMILES('[CH]=C=C([O])OC=C'),
    E0 = (249.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902698,0.0650968,-7.55594e-05,4.45932e-08,-1.02511e-11,30085.4,24.8396], Tmin=(100,'K'), Tmax=(1071.88,'K')), NASAPolynomial(coeffs=[14.3571,0.0148881,-5.29661e-06,8.92293e-10,-5.8485e-14,27201.1,-41.0038], Tmin=(1071.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CC(O)O[C]=C(23841)',
    structure = SMILES('[C]#CC(O)O[C]=C'),
    E0 = (459.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,347.232,347.232,347.232,347.232],'cm^-1')),
        HinderedRotor(inertia=(0.187029,'amu*angstrom^2'), symmetry=1, barrier=(16.002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.827904,'amu*angstrom^2'), symmetry=1, barrier=(70.8347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187029,'amu*angstrom^2'), symmetry=1, barrier=(16.002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187029,'amu*angstrom^2'), symmetry=1, barrier=(16.002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.553828,0.0730179,-9.71635e-05,6.37559e-08,-1.5909e-11,55447.8,27.7467], Tmin=(100,'K'), Tmax=(1066.64,'K')), NASAPolynomial(coeffs=[15.8488,0.0107077,-2.57285e-06,2.82144e-10,-1.17401e-14,52466.7,-45.7088], Tmin=(1066.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC(O)C#C(23842)',
    structure = SMILES('[CH]=[C]OC(O)C#C'),
    E0 = (369.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,1380,1390,370,380,2900,435,1685,370,2175,525,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.970452,'amu*angstrom^2'), symmetry=1, barrier=(22.3126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970361,'amu*angstrom^2'), symmetry=1, barrier=(22.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970399,'amu*angstrom^2'), symmetry=1, barrier=(22.3114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970403,'amu*angstrom^2'), symmetry=1, barrier=(22.3115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.581878,0.0714401,-8.5372e-05,4.5977e-08,-8.18851e-12,44617.7,27.3128], Tmin=(100,'K'), Tmax=(890.983,'K')), NASAPolynomial(coeffs=[17.459,0.00883257,-2.12707e-06,2.6864e-10,-1.5013e-14,41087.9,-55.0928], Tmin=(890.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[C]#CC([O])OC=C(23843)',
    structure = SMILES('[C]#CC([O])OC=C'),
    E0 = (445.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,306.239,306.585,306.766,306.983,307.679],'cm^-1')),
        HinderedRotor(inertia=(0.303975,'amu*angstrom^2'), symmetry=1, barrier=(20.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457036,'amu*angstrom^2'), symmetry=1, barrier=(30.5267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08958,'amu*angstrom^2'), symmetry=1, barrier=(72.6745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81917,0.0670877,-7.71118e-05,4.21966e-08,-8.05105e-12,53750,24.5323], Tmin=(100,'K'), Tmax=(892.604,'K')), NASAPolynomial(coeffs=[15.1242,0.0133113,-4.09806e-06,6.2726e-10,-3.89974e-14,50784.8,-45.1605], Tmin=(892.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=O(601)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,672.051,672.102],'cm^-1')),
        HinderedRotor(inertia=(0.000373196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57974,0.00389613,2.17609e-05,-3.06386e-08,1.18311e-11,19367.5,10.1348], Tmin=(100,'K'), Tmax=(961.532,'K')), NASAPolynomial(coeffs=[6.4326,0.00553733,-1.87382e-06,3.59985e-10,-2.76653e-14,18194.3,-6.76404], Tmin=(961.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][C]OC=O(4560)',
    structure = SMILES('[CH2][C]OC=O'),
    E0 = (264.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,243.595,243.598,1046.13],'cm^-1')),
        HinderedRotor(inertia=(0.0130986,'amu*angstrom^2'), symmetry=1, barrier=(25.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609379,'amu*angstrom^2'), symmetry=1, barrier=(25.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609381,'amu*angstrom^2'), symmetry=1, barrier=(25.6604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25739,0.0358204,-3.10362e-05,1.30264e-08,-2.17636e-12,31931.9,17.8429], Tmin=(100,'K'), Tmax=(1421.13,'K')), NASAPolynomial(coeffs=[10.818,0.0117254,-5.60416e-06,1.09602e-09,-7.76245e-14,29498.8,-26.4657], Tmin=(1421.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C=C([O])O[C]=C(23844)',
    structure = SMILES('[CH]=C=C([O])O[C]=C'),
    E0 = (488.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,1685,370,249.854,249.855,249.855],'cm^-1')),
        HinderedRotor(inertia=(0.226838,'amu*angstrom^2'), symmetry=1, barrier=(10.0489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75348,'amu*angstrom^2'), symmetry=1, barrier=(77.679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17517,0.0660994,-9.97567e-05,7.99103e-08,-2.50362e-11,58903.8,26.9917], Tmin=(100,'K'), Tmax=(872.632,'K')), NASAPolynomial(coeffs=[9.88621,0.0195824,-8.47406e-06,1.52258e-09,-1.00757e-13,57634.2,-12.4102], Tmin=(872.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]OC([O])C#C(23845)',
    structure = SMILES('[CH]=[C]OC([O])C#C'),
    E0 = (595.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,307.597,308.125,308.578],'cm^-1')),
        HinderedRotor(inertia=(0.265622,'amu*angstrom^2'), symmetry=1, barrier=(17.8372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380632,'amu*angstrom^2'), symmetry=1, barrier=(25.389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06011,'amu*angstrom^2'), symmetry=1, barrier=(70.9487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.97941,0.0682438,-9.59346e-05,6.85824e-08,-1.91281e-11,71744.2,26.7487], Tmin=(100,'K'), Tmax=(884.881,'K')), NASAPolynomial(coeffs=[12.6611,0.0154391,-6.42477e-06,1.14731e-09,-7.64523e-14,69676.8,-28.1804], Tmin=(884.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=CJO) + radical(CCOJ)"""),
)

species(
    label = '[C]#CC([O])O[C]=C(23846)',
    structure = SMILES('[C]#CC([O])O[C]=C'),
    E0 = (685.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,1685,370,278.992,280.2,280.417,280.499,1965.14],'cm^-1')),
        HinderedRotor(inertia=(1.41661,'amu*angstrom^2'), symmetry=1, barrier=(78.2683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285651,'amu*angstrom^2'), symmetry=1, barrier=(15.903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4136,'amu*angstrom^2'), symmetry=1, barrier=(78.2217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02118,0.0689561,-0.00010448,8.17882e-08,-2.47125e-11,82571.3,26.9348], Tmin=(100,'K'), Tmax=(913.938,'K')), NASAPolynomial(coeffs=[10.9901,0.0174245,-6.9368e-06,1.1769e-09,-7.45396e-14,81079.1,-18.4573], Tmin=(913.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(C=CJO) + radical(Acetyl)"""),
)

species(
    label = 'C=[C]OC1[C]=CO1(23807)',
    structure = SMILES('C=[C]OC1[C]=CO1'),
    E0 = (346.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975734,0.0440451,2.11243e-05,-8.08895e-08,4.01722e-11,41809.7,22.7835], Tmin=(100,'K'), Tmax=(923.668,'K')), NASAPolynomial(coeffs=[25.5975,-0.00431342,5.03261e-06,-9.79326e-10,5.87701e-14,34775.7,-107.502], Tmin=(923.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=C1C=[C]C([O])O1(23847)',
    structure = SMILES('C=C1C=[C]C([O])O1'),
    E0 = (212.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.952,0.028733,3.37124e-05,-7.00797e-08,3.02694e-11,25637.6,18.9327], Tmin=(100,'K'), Tmax=(952.991,'K')), NASAPolynomial(coeffs=[15.7964,0.0113922,-3.16278e-06,6.06212e-10,-4.98135e-14,21147.6,-56.905], Tmin=(952.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(=O)OC=C(23848)',
    structure = SMILES('C#CC(=O)OC=C'),
    E0 = (5.31192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07628,0.0658112,-8.60929e-05,5.85485e-08,-1.568e-11,742.841,22.2921], Tmin=(100,'K'), Tmax=(917.504,'K')), NASAPolynomial(coeffs=[12.1245,0.0176454,-7.34936e-06,1.33363e-09,-9.04106e-14,-1284.55,-30.0582], Tmin=(917.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.31192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1OC(=C)O1(23849)',
    structure = SMILES('C#CC1OC(=C)O1'),
    E0 = (72.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21251,0.0404106,2.55904e-05,-8.02399e-08,3.87216e-11,8856.8,16.4506], Tmin=(100,'K'), Tmax=(923.935,'K')), NASAPolynomial(coeffs=[22.9991,0.000181814,3.08363e-06,-6.34989e-10,3.64925e-14,2522.14,-99.4276], Tmin=(923.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
)

species(
    label = '[C]=C(584)',
    structure = SMILES('[C]=C'),
    E0 = (600.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94093,-0.00117598,1.80376e-05,-2.01208e-08,6.96659e-12,72197.9,5.25681], Tmin=(100,'K'), Tmax=(976.125,'K')), NASAPolynomial(coeffs=[3.93016,0.00536132,-1.98619e-06,3.69549e-10,-2.66221e-14,71890.7,3.724], Tmin=(976.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C#CC([O])[O](22347)',
    structure = SMILES('C#CC([O])[O]'),
    E0 = (263.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,957.767,2060.5],'cm^-1')),
        HinderedRotor(inertia=(2.63237,'amu*angstrom^2'), symmetry=1, barrier=(60.5233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63929,0.0399993,-8.12464e-05,8.74937e-08,-3.38057e-11,31722.1,16.5565], Tmin=(100,'K'), Tmax=(876.835,'K')), NASAPolynomial(coeffs=[0.362017,0.0251083,-1.25266e-05,2.36511e-09,-1.59456e-13,33093.3,32.7852], Tmin=(876.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=CO[C]=C(19367)',
    structure = SMILES('[CH]=C=CO[C]=C'),
    E0 = (504.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21977,'amu*angstrom^2'), symmetry=1, barrier=(28.0448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21614,'amu*angstrom^2'), symmetry=1, barrier=(27.9614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34459,0.0510081,-5.26791e-05,2.77943e-08,-5.63687e-12,60792.9,25.1839], Tmin=(100,'K'), Tmax=(1318.63,'K')), NASAPolynomial(coeffs=[13.3316,0.0108536,-2.68749e-06,3.38623e-10,-1.80043e-14,57961.3,-34.7116], Tmin=(1318.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (348.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (432.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (445.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (489.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (545.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (348.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (555.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (502.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (461.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (490.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (543.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (402.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (490.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (430.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (816.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (700.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (807.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (897.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (590.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (407.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (437.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (356.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (868.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (911.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C=C=O(598)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['[CH]=C1OC1O[C]=C(23834)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['[CH]=C1C(=C)OC1[O](23835)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.06165e+08,'s^-1'), n=1.17712, Ea=(96.5377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(=O)O[C]=C(23836)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#COC([O])C#C(23837)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C=O(598)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(138.052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from 134.2 to 138.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]#C(5143)', 'C=[C]OC=O(4555)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-NdH_O;CtJ_Ct]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C#C[C](O)O[C]=C(23838)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_Ct]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=COC([O])C#C(23839)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['[CH]=C=C([O])OC=C(23840)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]#CC(O)O[C]=C(23841)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC(O)C#C(23842)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC([O])OC=C(23843)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(601)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#C(5143)', '[CH2][C]OC=O(4560)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]=C=C([O])O[C]=C(23844)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]=[C]OC([O])C#C(23845)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[C]#CC([O])O[C]=C(23846)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C=[C]OC1[C]=CO1(23807)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C=C1C=[C]C([O])O1(23847)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for R5_SS_T;triplebond_intra_H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C#CC(=O)OC=C(23848)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([O])O[C]=C(22460)'],
    products = ['C#CC1OC(=C)O1(23849)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=C(584)', 'C#CC([O])[O](22347)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(T)(63)', '[CH]=C=CO[C]=C(19367)'],
    products = ['C#CC([O])O[C]=C(22460)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

network(
    label = '4757',
    isomers = [
        'C#CC([O])O[C]=C(22460)',
    ],
    reactants = [
        ('C=C=O(598)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4757',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)


species(
    label = 'C#CC([O])CC=C[O](22586)',
    structure = SMILES('C#CC([O])CC=C[O]'),
    E0 = (219.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,218.439,218.439,218.439,218.439],'cm^-1')),
        HinderedRotor(inertia=(0.801212,'amu*angstrom^2'), symmetry=1, barrier=(27.1292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801212,'amu*angstrom^2'), symmetry=1, barrier=(27.1292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.801212,'amu*angstrom^2'), symmetry=1, barrier=(27.1292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349841,0.0658846,-3.00643e-05,-2.26308e-08,1.7727e-11,26589.5,30.1039], Tmin=(100,'K'), Tmax=(927.933,'K')), NASAPolynomial(coeffs=[21.5129,0.0116462,-2.17941e-06,2.92247e-10,-2.20858e-14,21069.4,-78.9933], Tmin=(927.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC=O(5269)',
    structure = SMILES('C=CC=O'),
    E0 = (-81.3387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.873408,'amu*angstrom^2'), symmetry=1, barrier=(20.0814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3136.31,'J/mol'), sigma=(5.14154,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.88 K, Pc=52.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.9738,0.0193269,-1.02836e-06,-7.40922e-09,2.6466e-12,-9743.32,12.1361], Tmin=(100,'K'), Tmax=(1315.19,'K')), NASAPolynomial(coeffs=[7.40832,0.0154746,-7.62321e-06,1.50372e-09,-1.06406e-13,-11743,-13.6408], Tmin=(1315.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.3387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
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
    label = '[CH]=C1OC1CC=C[O](25778)',
    structure = SMILES('[CH]=C1OC1CC=C[O]'),
    E0 = (241.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361482,0.0536966,2.62341e-05,-9.63734e-08,4.74296e-11,29256.1,26.8998], Tmin=(100,'K'), Tmax=(924.561,'K')), NASAPolynomial(coeffs=[28.7499,-0.000973057,4.36488e-06,-8.80087e-10,5.10773e-14,21094,-123.583], Tmin=(924.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C1CC(C=O)O1(25421)',
    structure = SMILES('[CH]=[C]C1CC(C=O)O1'),
    E0 = (339.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40459,0.0399982,3.1727e-05,-7.9595e-08,3.70866e-11,40909,28.8871], Tmin=(100,'K'), Tmax=(907.328,'K')), NASAPolynomial(coeffs=[17.2426,0.0153423,-2.1813e-06,1.83403e-10,-1.19745e-14,36175.8,-56.2269], Tmin=(907.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC=CCC1[O](25509)',
    structure = SMILES('[CH]=C1OC=CCC1[O]'),
    E0 = (248.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49029,0.030555,6.96985e-05,-1.14603e-07,4.56632e-11,30049.8,23.5096], Tmin=(100,'K'), Tmax=(981.579,'K')), NASAPolynomial(coeffs=[19.0646,0.0201549,-7.95621e-06,1.67401e-09,-1.33553e-13,23650.6,-75.9715], Tmin=(981.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(CC(C)OJ) + radical(Cds_P)"""),
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
    label = 'C#CC(=O)CC=C[O](25779)',
    structure = SMILES('C#CC(=O)CC=C[O]'),
    E0 = (56.5766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,375,552.5,462.5,1710,277.611,277.611,277.611],'cm^-1')),
        HinderedRotor(inertia=(0.498994,'amu*angstrom^2'), symmetry=1, barrier=(27.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498993,'amu*angstrom^2'), symmetry=1, barrier=(27.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498994,'amu*angstrom^2'), symmetry=1, barrier=(27.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69711,0.061702,-4.2e-05,7.28591e-09,1.83097e-12,6932.75,25.898], Tmin=(100,'K'), Tmax=(1117.15,'K')), NASAPolynomial(coeffs=[17.5688,0.019125,-8.77592e-06,1.74822e-09,-1.27456e-13,2050.33,-62.3477], Tmin=(1117.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.5766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])CC=C=O(25780)',
    structure = SMILES('C#CC([O])CC=C=O'),
    E0 = (202.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2175,525,1380,1390,370,380,2900,435,345.718,345.718,345.722],'cm^-1')),
        HinderedRotor(inertia=(0.22052,'amu*angstrom^2'), symmetry=1, barrier=(18.7049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220546,'amu*angstrom^2'), symmetry=1, barrier=(18.7049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35332,'amu*angstrom^2'), symmetry=1, barrier=(114.794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627362,0.0696417,-7.14791e-05,3.78738e-08,-7.9253e-12,24425.4,28.5742], Tmin=(100,'K'), Tmax=(1165.8,'K')), NASAPolynomial(coeffs=[15.1001,0.0199841,-7.58643e-06,1.33666e-09,-9.01275e-14,21050.9,-43.4686], Tmin=(1165.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C[O](5266)',
    structure = SMILES('[CH2]C=C[O]'),
    E0 = (90.2929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.57685,'amu*angstrom^2'), symmetry=1, barrier=(36.2549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69019,0.0144913,4.15491e-05,-7.27602e-08,3.14101e-11,10920.2,13.4175], Tmin=(100,'K'), Tmax=(922.751,'K')), NASAPolynomial(coeffs=[14.044,0.00224417,1.35973e-06,-3.04875e-10,1.62832e-14,7250.86,-48.974], Tmin=(922.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.2929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CCC=O(4687)',
    structure = SMILES('[O]C=CCC=O'),
    E0 = (-160.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,242.784,243.031],'cm^-1')),
        HinderedRotor(inertia=(0.547197,'amu*angstrom^2'), symmetry=1, barrier=(23.8946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545022,'amu*angstrom^2'), symmetry=1, barrier=(23.8969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94847,0.0299836,2.45267e-05,-5.53063e-08,2.30245e-11,-19212.7,20.3002], Tmin=(100,'K'), Tmax=(1003.7,'K')), NASAPolynomial(coeffs=[15.111,0.014027,-6.17416e-06,1.31662e-09,-1.03456e-13,-23693.4,-52.4087], Tmin=(1003.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ)"""),
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
    label = '[CH]=C=C(O)CC=C[O](25781)',
    structure = SMILES('[CH]=C=C(O)CC=C[O]'),
    E0 = (101.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36492,0.0877476,-9.17679e-05,4.48091e-08,-7.9882e-12,12477,34.9153], Tmin=(100,'K'), Tmax=(1643.99,'K')), NASAPolynomial(coeffs=[24.3601,0.00510607,1.92915e-06,-6.05141e-10,4.60302e-14,6728.06,-93.7414], Tmin=(1643.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([O])CC=[C]O(25782)',
    structure = SMILES('C#CC([O])CC=[C]O'),
    E0 = (318.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2175,525,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,260.068,260.081,260.089],'cm^-1')),
        HinderedRotor(inertia=(0.434159,'amu*angstrom^2'), symmetry=1, barrier=(20.8393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434166,'amu*angstrom^2'), symmetry=1, barrier=(20.8393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434174,'amu*angstrom^2'), symmetry=1, barrier=(20.8394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434137,'amu*angstrom^2'), symmetry=1, barrier=(20.8392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.515749,0.0824061,-8.74718e-05,4.49883e-08,-8.68006e-12,38442.4,34.9308], Tmin=(100,'K'), Tmax=(1440.43,'K')), NASAPolynomial(coeffs=[21.7676,0.0100353,-1.18295e-06,-4.71542e-12,6.48885e-15,33111.2,-76.9276], Tmin=(1440.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(O)[CH]C=C[O](25783)',
    structure = SMILES('C#CC(O)[CH]C=C[O]'),
    E0 = (106.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0238042,0.0737386,-4.33502e-05,-1.4298e-08,1.59966e-11,12959.1,28.8364], Tmin=(100,'K'), Tmax=(924.681,'K')), NASAPolynomial(coeffs=[23.6945,0.00944229,-1.18738e-06,1.03077e-10,-8.96117e-15,6935.14,-92.5885], Tmin=(924.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C#CC([O])C[C]=CO(25784)',
    structure = SMILES('C#CC([O])C[C]=CO'),
    E0 = (316.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2175,525,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04925,'amu*angstrom^2'), symmetry=1, barrier=(24.1242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04952,'amu*angstrom^2'), symmetry=1, barrier=(24.1306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04934,'amu*angstrom^2'), symmetry=1, barrier=(24.1264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04977,'amu*angstrom^2'), symmetry=1, barrier=(24.1363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02082,0.0883127,-9.60176e-05,4.92536e-08,-9.32752e-12,38236.2,34.0908], Tmin=(100,'K'), Tmax=(1507.81,'K')), NASAPolynomial(coeffs=[24.1261,0.00618192,1.02745e-06,-4.36437e-10,3.58037e-14,32405.7,-91.7432], Tmin=(1507.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(O)C[C]=C[O](25785)',
    structure = SMILES('C#CC(O)C[C]=C[O]'),
    E0 = (227.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2175,525,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,260.583,260.824,260.929],'cm^-1')),
        HinderedRotor(inertia=(0.447391,'amu*angstrom^2'), symmetry=1, barrier=(21.6246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44804,'amu*angstrom^2'), symmetry=1, barrier=(21.6239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447671,'amu*angstrom^2'), symmetry=1, barrier=(21.6225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.448145,'amu*angstrom^2'), symmetry=1, barrier=(21.623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.458364,0.0816396,-8.55559e-05,4.35494e-08,-8.35103e-12,27519.1,33.6736], Tmin=(100,'K'), Tmax=(1435.42,'K')), NASAPolynomial(coeffs=[21.6664,0.0107077,-1.73722e-06,1.17667e-10,-2.45997e-15,22123.2,-77.7335], Tmin=(1435.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])[CH]C=CO(25786)',
    structure = SMILES('C#CC([O])[CH]C=CO'),
    E0 = (195.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75605,0.0938382,-9.92326e-05,4.84379e-08,-8.61675e-12,23728.1,33.4771], Tmin=(100,'K'), Tmax=(1646.74,'K')), NASAPolynomial(coeffs=[26.4132,0.00381221,2.44817e-06,-6.92597e-10,5.13729e-14,17379.5,-107.581], Tmin=(1646.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC(O)CC=C[O](25787)',
    structure = SMILES('[C]#CC(O)CC=C[O]'),
    E0 = (326.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,335.926,336.024,336.418,336.938],'cm^-1')),
        HinderedRotor(inertia=(0.217919,'amu*angstrom^2'), symmetry=1, barrier=(17.6772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220576,'amu*angstrom^2'), symmetry=1, barrier=(17.7078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218091,'amu*angstrom^2'), symmetry=1, barrier=(17.6923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15195,'amu*angstrom^2'), symmetry=1, barrier=(92.6287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437916,0.0810476,-8.51962e-05,4.3912e-08,-8.50269e-12,39461.7,33.8656], Tmin=(100,'K'), Tmax=(1441.42,'K')), NASAPolynomial(coeffs=[20.8216,0.0113188,-1.46471e-06,1.98293e-11,5.91423e-15,34447.9,-72.605], Tmin=(1441.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=COJ)"""),
)

species(
    label = 'C#CC(O)C[CH][C]=O(25788)',
    structure = SMILES('C#CC(O)C[CH][C]=O'),
    E0 = (173.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299969,0.077542,-8.77427e-05,5.1858e-08,-1.20064e-11,21034,31.3407], Tmin=(100,'K'), Tmax=(1063.08,'K')), NASAPolynomial(coeffs=[15.5465,0.0201757,-6.8005e-06,1.09925e-09,-6.98701e-14,17792.3,-43.1476], Tmin=(1063.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = '[CH]=C=C([O])CC=CO(25789)',
    structure = SMILES('[CH]=C=C([O])CC=CO'),
    E0 = (98.2588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2211,0.0887642,-9.68632e-05,4.9591e-08,-9.27175e-12,12028.2,34.4727], Tmin=(100,'K'), Tmax=(1561.89,'K')), NASAPolynomial(coeffs=[23.9975,0.00477969,2.42451e-06,-7.40745e-10,5.73532e-14,6516.78,-90.862], Tmin=(1561.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.2588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])CC=CO(25790)',
    structure = SMILES('[C]#CC([O])CC=CO'),
    E0 = (415.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,197.294,197.42,197.545,197.768],'cm^-1')),
        HinderedRotor(inertia=(0.81616,'amu*angstrom^2'), symmetry=1, barrier=(22.6351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819593,'amu*angstrom^2'), symmetry=1, barrier=(22.6332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81616,'amu*angstrom^2'), symmetry=1, barrier=(22.631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817339,'amu*angstrom^2'), symmetry=1, barrier=(22.6331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.998345,0.0877007,-9.56042e-05,4.95654e-08,-9.46396e-12,50178.8,34.2752], Tmin=(100,'K'), Tmax=(1509.25,'K')), NASAPolynomial(coeffs=[23.2098,0.00689682,1.24677e-06,-5.2279e-10,4.32893e-14,44767.2,-86.1994], Tmin=(1509.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl)"""),
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
    label = '[CH]=C[O](602)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (221.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27415,0.00479611,3.742e-05,-6.11894e-08,2.65903e-11,26726.9,9.63858], Tmin=(100,'K'), Tmax=(905.806,'K')), NASAPolynomial(coeffs=[11.9892,-0.00434473,3.96329e-06,-8.00891e-10,5.23184e-14,23944.2,-38.1893], Tmin=(905.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([CH2])[O](5294)',
    structure = SMILES('C#CC([CH2])[O]'),
    E0 = (418.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3100,440,815,1455,1000,348.129],'cm^-1')),
        HinderedRotor(inertia=(0.289811,'amu*angstrom^2'), symmetry=1, barrier=(24.9417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397394,'amu*angstrom^2'), symmetry=1, barrier=(34.187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46647,0.0444261,-4.38277e-05,2.13721e-08,-3.91288e-12,50443.4,20.5021], Tmin=(100,'K'), Tmax=(1545.21,'K')), NASAPolynomial(coeffs=[13.2214,0.00716331,-1.02144e-06,4.17503e-11,1.21996e-15,47626.5,-38.6839], Tmin=(1545.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C#CC([O])[CH]C=C[O](25791)',
    structure = SMILES('C#CC([O])[CH]C=C[O]'),
    E0 = (336.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.70218,'amu*angstrom^2'), symmetry=1, barrier=(39.1364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70131,'amu*angstrom^2'), symmetry=1, barrier=(39.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70144,'amu*angstrom^2'), symmetry=1, barrier=(39.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.13379,0.0698092,-3.76847e-05,-1.89953e-08,1.74448e-11,40659.9,28.0787], Tmin=(100,'K'), Tmax=(926.931,'K')), NASAPolynomial(coeffs=[23.6379,0.00763501,-5.93267e-07,1.36443e-11,-3.81507e-15,34616.2,-92.6278], Tmin=(926.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O][CH]CC=C[O](1115)',
    structure = SMILES('[O][CH]CC=C[O]'),
    E0 = (162.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,180,180,525.062],'cm^-1')),
        HinderedRotor(inertia=(0.0297468,'amu*angstrom^2'), symmetry=1, barrier=(12.4467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0421597,'amu*angstrom^2'), symmetry=1, barrier=(16.9616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56518,0.0527261,-5.16335e-05,2.68534e-08,-5.65794e-12,19611.4,23.2023], Tmin=(100,'K'), Tmax=(1139.15,'K')), NASAPolynomial(coeffs=[10.9093,0.019915,-8.42836e-06,1.5682e-09,-1.0876e-13,17482.6,-23.095], Tmin=(1139.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCsJOH) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=C([O])CC=C[O](25792)',
    structure = SMILES('[CH]=C=C([O])CC=C[O]'),
    E0 = (239.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03126,'amu*angstrom^2'), symmetry=1, barrier=(23.7107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0317,'amu*angstrom^2'), symmetry=1, barrier=(23.7208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.553846,0.0787383,-8.25523e-05,4.12837e-08,-7.63531e-12,29014.2,33.4906], Tmin=(100,'K'), Tmax=(1550.66,'K')), NASAPolynomial(coeffs=[21.6072,0.00729702,3.64395e-07,-3.01194e-10,2.62585e-14,23857.7,-77.6112], Tmin=(1550.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([O])C[C]=C[O](25793)',
    structure = SMILES('C#CC([O])C[C]=C[O]'),
    E0 = (457.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,2175,525,241.484,241.507,241.511,241.533],'cm^-1')),
        HinderedRotor(inertia=(0.647163,'amu*angstrom^2'), symmetry=1, barrier=(26.7827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647104,'amu*angstrom^2'), symmetry=1, barrier=(26.7828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647005,'amu*angstrom^2'), symmetry=1, barrier=(26.7829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.353346,0.0782859,-8.17103e-05,4.09578e-08,-7.69754e-12,55222.2,33.1077], Tmin=(100,'K'), Tmax=(1477.61,'K')), NASAPolynomial(coeffs=[21.5698,0.00893296,-1.14951e-06,2.78721e-11,2.82154e-15,49835.7,-77.5222], Tmin=(1477.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])C[CH][C]=O(25794)',
    structure = SMILES('C#CC([O])C[CH][C]=O'),
    E0 = (404.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2175,525,418.555,419.413,419.998],'cm^-1')),
        HinderedRotor(inertia=(0.0830488,'amu*angstrom^2'), symmetry=1, barrier=(10.361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000959974,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254443,'amu*angstrom^2'), symmetry=1, barrier=(31.5658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669731,'amu*angstrom^2'), symmetry=1, barrier=(83.3619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444524,0.0737631,-8.25833e-05,4.77777e-08,-1.08e-11,48735.4,30.63], Tmin=(100,'K'), Tmax=(1089.37,'K')), NASAPolynomial(coeffs=[15.5648,0.0182433,-6.1353e-06,9.93197e-10,-6.33568e-14,45441.1,-43.6109], Tmin=(1089.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCHO) + radical(CCCJ=O) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])CC=C[O](25795)',
    structure = SMILES('[C]#CC([O])CC=C[O]'),
    E0 = (557.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,332.657,332.663,332.666,332.67,332.685],'cm^-1')),
        HinderedRotor(inertia=(0.284382,'amu*angstrom^2'), symmetry=1, barrier=(22.3286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97049,'amu*angstrom^2'), symmetry=1, barrier=(76.2166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284321,'amu*angstrom^2'), symmetry=1, barrier=(22.329,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331805,0.0776831,-8.13217e-05,4.12929e-08,-7.84102e-12,67164.8,33.2956], Tmin=(100,'K'), Tmax=(1481.14,'K')), NASAPolynomial(coeffs=[20.6842,0.00960342,-9.07492e-07,-6.33715e-11,1.06849e-14,62181.3,-72.1568], Tmin=(1481.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])CC1[CH]O1(25796)',
    structure = SMILES('C#CC([O])CC1[CH]O1'),
    E0 = (360.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.784952,0.0830335,-8.76255e-05,4.59597e-08,-8.86604e-12,43518,32.0277], Tmin=(100,'K'), Tmax=(1529.4,'K')), NASAPolynomial(coeffs=[19.23,0.011879,6.07158e-07,-5.41653e-10,4.95503e-14,39595.4,-65.8456], Tmin=(1529.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=CCC1[C]=CO1(25704)',
    structure = SMILES('[O]C=CCC1[C]=CO1'),
    E0 = (213.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5358,0.0456415,5.52791e-05,-1.29128e-07,5.95499e-11,25799.5,26.9175], Tmin=(100,'K'), Tmax=(925.362,'K')), NASAPolynomial(coeffs=[29.9136,-0.00246227,5.38206e-06,-1.05575e-09,6.06286e-14,16985,-130.784], Tmin=(925.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC1C[CH]C([O])O1(25797)',
    structure = SMILES('C#CC1C[CH]C([O])O1'),
    E0 = (260.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7312,0.0347579,3.8881e-05,-8.18284e-08,3.70618e-11,31367.5,24.6301], Tmin=(100,'K'), Tmax=(895.433,'K')), NASAPolynomial(coeffs=[14.2348,0.0194913,-3.53688e-06,3.73793e-10,-2.18627e-14,27501.1,-43.3974], Tmin=(895.433,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[C]=COC=CC1(25563)',
    structure = SMILES('[O]C1[C]=COC=CC1'),
    E0 = (275.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39263,-0.00576766,0.000250293,-3.74275e-07,1.59471e-10,33299.7,27.4643], Tmin=(100,'K'), Tmax=(900.618,'K')), NASAPolynomial(coeffs=[46.1745,-0.0387217,2.88021e-05,-5.73722e-09,3.79699e-13,18503.6,-221.257], Tmin=(900.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(=O)CC=CO(25798)',
    structure = SMILES('C#CC(=O)CC=CO'),
    E0 = (-84.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.304075,0.0685628,-4.55773e-05,2.18623e-09,5.67212e-12,-10065.2,25.892], Tmin=(100,'K'), Tmax=(1022.41,'K')), NASAPolynomial(coeffs=[20.1343,0.0165444,-6.76512e-06,1.33385e-09,-9.92363e-14,-15456.3,-76.751], Tmin=(1022.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(O)CC=C=O(25799)',
    structure = SMILES('C#CC(O)CC=C=O'),
    E0 = (-28.3183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493729,0.0732941,-7.62116e-05,4.14301e-08,-8.92342e-12,-3276.37,29.2457], Tmin=(100,'K'), Tmax=(1132.25,'K')), NASAPolynomial(coeffs=[14.9974,0.022055,-8.32959e-06,1.46079e-09,-9.81205e-14,-6560.69,-42.5274], Tmin=(1132.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.3183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([O])[CH]C[CH][O](25800)',
    structure = SMILES('C#CC([O])[CH]C[CH][O]'),
    E0 = (605.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404706,0.0832451,-0.00011084,8.22945e-08,-2.46703e-11,72950.4,33.7797], Tmin=(100,'K'), Tmax=(814.571,'K')), NASAPolynomial(coeffs=[11.1933,0.0302679,-1.32861e-05,2.45553e-09,-1.67406e-13,71192.7,-16.0568], Tmin=(814.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C=C([O])CC[CH][O](25801)',
    structure = SMILES('[CH]=C=C([O])CC[CH][O]'),
    E0 = (424.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,350,440,435,1725,3025,407.5,1350,352.5,267.941,267.949,267.951,3289.73],'cm^-1')),
        HinderedRotor(inertia=(0.164077,'amu*angstrom^2'), symmetry=1, barrier=(8.35722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164008,'amu*angstrom^2'), symmetry=1, barrier=(8.3572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47581,'amu*angstrom^2'), symmetry=1, barrier=(75.191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125977,0.0919654,-0.000140233,1.15684e-07,-3.71887e-11,51149.6,32.3864], Tmin=(100,'K'), Tmax=(881.843,'K')), NASAPolynomial(coeffs=[10.6464,0.0307798,-1.32536e-05,2.37265e-09,-1.56397e-13,49817.6,-14.0772], Tmin=(881.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(CCsJOH) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=C([O])C[CH]C[O](25802)',
    structure = SMILES('[CH]=C=C([O])C[CH]C[O]'),
    E0 = (443.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,350,440,435,1725,3025,407.5,1350,352.5,257.23,276.583,282.701,1838.31],'cm^-1')),
        HinderedRotor(inertia=(0.00240543,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148182,'amu*angstrom^2'), symmetry=1, barrier=(7.6214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215913,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611559,0.0777047,-9.90824e-05,7.18321e-08,-2.11049e-11,53493.4,33.5944], Tmin=(100,'K'), Tmax=(830.326,'K')), NASAPolynomial(coeffs=[10.5706,0.0297268,-1.24067e-05,2.23828e-09,-1.50499e-13,51839.6,-12.6], Tmin=(830.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[C]#CC([O])CC[CH][O](25803)',
    structure = SMILES('[C]#CC([O])CC[CH][O]'),
    E0 = (742.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2175,525,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145571,0.0931172,-0.000146334,1.24794e-07,-4.11637e-11,89462,32.9498], Tmin=(100,'K'), Tmax=(885.27,'K')), NASAPolynomial(coeffs=[9.37392,0.0334056,-1.46354e-05,2.63035e-09,-1.73346e-13,88534,-6.46017], Tmin=(885.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[C]#CC([O])C[CH]C[O](25804)',
    structure = SMILES('[C]#CC([O])C[CH]C[O]'),
    E0 = (762.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2175,525,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595724,0.0793126,-0.00010697,8.35595e-08,-2.63435e-11,91807.3,34.2825], Tmin=(100,'K'), Tmax=(850.506,'K')), NASAPolynomial(coeffs=[9.3894,0.0321892,-1.3691e-05,2.47232e-09,-1.65448e-13,90520,-5.49202], Tmin=(850.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCO) + radical(CC(C)OJ) + radical(CCOJ) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1CC=COO1(25617)',
    structure = SMILES('C#CC1CC=COO1'),
    E0 = (195.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22818,0.0415511,3.57251e-05,-8.86586e-08,4.13433e-11,23620.5,21.4418], Tmin=(100,'K'), Tmax=(909.647,'K')), NASAPolynomial(coeffs=[19.455,0.0125558,-8.13909e-07,-5.93709e-11,3.21474e-15,18188.2,-76.3989], Tmin=(909.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(34dihydro12dioxin)"""),
)

species(
    label = 'C=CC([O])C=C=C[O](22584)',
    structure = SMILES('C=CC([O])C=C=C[O]'),
    E0 = (214.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438181,0.0654159,-3.65856e-05,-7.54109e-09,9.51967e-12,25941.3,30.6225], Tmin=(100,'K'), Tmax=(987.953,'K')), NASAPolynomial(coeffs=[19.5037,0.016784,-6.11062e-06,1.15527e-09,-8.53236e-14,20780.3,-68.1804], Tmin=(987.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[CH]CC([O])[C]=C1(25805)',
    structure = SMILES('[O]C1[CH]CC([O])[C]=C1'),
    E0 = (531.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47703,0.0342539,4.97338e-05,-9.31754e-08,3.86857e-11,64044.6,27.9993], Tmin=(100,'K'), Tmax=(969.42,'K')), NASAPolynomial(coeffs=[18.2209,0.0177717,-6.16096e-06,1.24023e-09,-9.85089e-14,58326.4,-65.0093], Tmin=(969.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(Cds_S) + radical(CC(C)OJ) + radical(CCJCO)"""),
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
    label = 'C#C[CH]CC=C[O](20022)',
    structure = SMILES('C#C[CH]CC=C[O]'),
    E0 = (308.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,3025,407.5,1350,352.5,287.209,289.294,289.85],'cm^-1')),
        HinderedRotor(inertia=(0.47549,'amu*angstrom^2'), symmetry=1, barrier=(27.955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59061,'amu*angstrom^2'), symmetry=1, barrier=(34.3175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12682,'amu*angstrom^2'), symmetry=1, barrier=(66.5005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0647921,0.0677002,-6.37812e-05,3.0016e-08,-5.2724e-12,37291.1,28.5977], Tmin=(100,'K'), Tmax=(1651.61,'K')), NASAPolynomial(coeffs=[17.2967,0.0123796,-1.48396e-06,4.11297e-12,6.9518e-15,33366.6,-58.3917], Tmin=(1651.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC([O])C#C(24341)',
    structure = SMILES('[CH]=CCC([O])C#C'),
    E0 = (534.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2175,525,317.184,317.496],'cm^-1')),
        HinderedRotor(inertia=(0.310077,'amu*angstrom^2'), symmetry=1, barrier=(22.236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309986,'amu*angstrom^2'), symmetry=1, barrier=(22.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310755,'amu*angstrom^2'), symmetry=1, barrier=(22.2201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883734,0.060839,-4.64472e-05,1.08371e-08,2.05497e-12,64380.6,27.3205], Tmin=(100,'K'), Tmax=(967.862,'K')), NASAPolynomial(coeffs=[15.3643,0.0178823,-6.04734e-06,1.03887e-09,-7.11148e-14,60786.5,-46.1536], Tmin=(967.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1C([O])CC1C=O(25339)',
    structure = SMILES('[CH]=C1C([O])CC1C=O'),
    E0 = (291.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.5053,-0.0122532,0.000119821,-1.16406e-07,2.65028e-11,34658.7,-18.2991], Tmin=(100,'K'), Tmax=(1757.72,'K')), NASAPolynomial(coeffs=[82.9578,0.0164676,-6.7961e-05,1.67419e-08,-1.24244e-12,-20015.4,-487.39], Tmin=(1757.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([O])C=CC=O(25806)',
    structure = SMILES('C#CC([O])C=CC=O'),
    E0 = (183.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.97645,'amu*angstrom^2'), symmetry=1, barrier=(22.4505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975533,'amu*angstrom^2'), symmetry=1, barrier=(22.4294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975797,'amu*angstrom^2'), symmetry=1, barrier=(22.4355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1004,0.065935,-6.36436e-05,3.21031e-08,-6.66085e-12,22201.3,26.8158], Tmin=(100,'K'), Tmax=(1138.87,'K')), NASAPolynomial(coeffs=[11.9924,0.0276795,-1.32574e-05,2.6082e-09,-1.86236e-13,19720.4,-27.1481], Tmin=(1138.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH]CC=O(25486)',
    structure = SMILES('C#CC([O])[CH]CC=O'),
    E0 = (276.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2175,525,180,180.183,908.822],'cm^-1')),
        HinderedRotor(inertia=(0.0028609,'amu*angstrom^2'), symmetry=1, barrier=(1.67803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638009,'amu*angstrom^2'), symmetry=1, barrier=(14.6691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638593,'amu*angstrom^2'), symmetry=1, barrier=(14.6827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133412,'amu*angstrom^2'), symmetry=1, barrier=(78.0183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723987,0.0698146,-6.86251e-05,3.5932e-08,-7.58446e-12,33378.1,31.547], Tmin=(100,'K'), Tmax=(1142.25,'K')), NASAPolynomial(coeffs=[13.4186,0.0253598,-1.0247e-05,1.85989e-09,-1.27199e-13,30478,-31.3854], Tmin=(1142.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])CC[C]=O(22575)',
    structure = SMILES('C#CC([O])CC[C]=O'),
    E0 = (236.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,316.028,316.556,317.245],'cm^-1')),
        HinderedRotor(inertia=(0.00167425,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163237,'amu*angstrom^2'), symmetry=1, barrier=(11.6267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16276,'amu*angstrom^2'), symmetry=1, barrier=(11.6311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.443249,'amu*angstrom^2'), symmetry=1, barrier=(31.5687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520579,0.0771869,-8.92701e-05,5.55555e-08,-1.38336e-11,28579.1,30.6033], Tmin=(100,'K'), Tmax=(978.946,'K')), NASAPolynomial(coeffs=[13.1071,0.0257573,-1.04658e-05,1.88872e-09,-1.28215e-13,26114.9,-29.8515], Tmin=(978.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C=C([O])CCC=O(25807)',
    structure = SMILES('[CH]=C=C([O])CCC=O'),
    E0 = (95.2002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583737,0.0767633,-9.12296e-05,5.97968e-08,-1.57753e-11,11571.5,29.6667], Tmin=(100,'K'), Tmax=(923.886,'K')), NASAPolynomial(coeffs=[11.9792,0.0274263,-1.11271e-05,1.99554e-09,-1.34501e-13,9465.93,-24.4074], Tmin=(923.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.2002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CC([O])CCC=O(25808)',
    structure = SMILES('[C]#CC([O])CCC=O'),
    E0 = (413.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2175,525,2782.5,750,1395,475,1775,1000,355.461,355.461,355.461,355.461],'cm^-1')),
        HinderedRotor(inertia=(0.00133418,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0804519,'amu*angstrom^2'), symmetry=1, barrier=(7.21354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0804521,'amu*angstrom^2'), symmetry=1, barrier=(7.21355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72585,'amu*angstrom^2'), symmetry=1, barrier=(65.0818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583171,0.0781568,-9.81994e-05,7.00882e-08,-2.02885e-11,49884.8,30.3022], Tmin=(100,'K'), Tmax=(842.202,'K')), NASAPolynomial(coeffs=[10.7129,0.0300444,-1.25059e-05,2.25276e-09,-1.51427e-13,48178.6,-16.8277], Tmin=(842.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[C]=CC(C=O)C1(25389)',
    structure = SMILES('[O]C1[C]=CC(C=O)C1'),
    E0 = (223.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32417,0.0476189,-7.64945e-06,-2.15865e-08,1.0607e-11,26996.8,25.9062], Tmin=(100,'K'), Tmax=(1053.27,'K')), NASAPolynomial(coeffs=[13.3816,0.025845,-1.08434e-05,2.08397e-09,-1.49807e-13,23124.7,-39.2136], Tmin=(1053.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C#CC1C[CH][CH]OO1(25601)',
    structure = SMILES('C#CC1C[CH][CH]OO1'),
    E0 = (463.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12276,0.0461306,2.52131e-05,-8.13817e-08,4.07386e-11,55877.3,24.685], Tmin=(100,'K'), Tmax=(879.413,'K')), NASAPolynomial(coeffs=[19.1335,0.0123689,6.54705e-07,-4.91662e-10,4.02156e-14,50847.2,-70.4799], Tmin=(879.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxane) + radical(CCsJOOC) + radical(CCJCOOH)"""),
)

species(
    label = 'C#CC(O)C=CC=O(25809)',
    structure = SMILES('C#CC(O)C=CC=O'),
    E0 = (-46.6214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.945021,0.0698098,-6.90141e-05,3.6327e-08,-7.88966e-12,-5499.41,27.5678], Tmin=(100,'K'), Tmax=(1091.12,'K')), NASAPolynomial(coeffs=[11.851,0.0298289,-1.40513e-05,2.74529e-09,-1.95364e-13,-7879.38,-25.9985], Tmin=(1091.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(=O)CCC=O(25810)',
    structure = SMILES('C#CC(=O)CCC=O'),
    E0 = (-92.8574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770846,0.0807783,-0.000126251,1.1734e-07,-4.28928e-11,-11061.3,26.3573], Tmin=(100,'K'), Tmax=(834.331,'K')), NASAPolynomial(coeffs=[4.70752,0.0418663,-2.02668e-05,3.86722e-09,-2.66095e-13,-11020.8,12.2578], Tmin=(834.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.8574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([O])C([CH2])C=O(22594)',
    structure = SMILES('C#CC([O])C([CH2])C=O'),
    E0 = (281.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,388.962,388.962],'cm^-1')),
        HinderedRotor(inertia=(0.0974447,'amu*angstrom^2'), symmetry=1, barrier=(10.4616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974439,'amu*angstrom^2'), symmetry=1, barrier=(10.4616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197874,'amu*angstrom^2'), symmetry=1, barrier=(21.2436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674186,'amu*angstrom^2'), symmetry=1, barrier=(72.3796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4180.37,'J/mol'), sigma=(6.7983,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=652.96 K, Pc=30.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337006,0.0819137,-9.87928e-05,6.33783e-08,-1.6199e-11,33972.1,30.6582], Tmin=(100,'K'), Tmax=(955.812,'K')), NASAPolynomial(coeffs=[13.7455,0.0257999,-1.07305e-05,1.95583e-09,-1.33398e-13,31408.9,-33.4241], Tmin=(955.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#CC1CC(C=O)O1(22600)',
    structure = SMILES('C#CC1CC(C=O)O1'),
    E0 = (20.2394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49021,0.0346665,5.47604e-05,-1.10589e-07,5.06067e-11,2544.02,25.3224], Tmin=(100,'K'), Tmax=(890.104,'K')), NASAPolynomial(coeffs=[19.1387,0.0108414,1.40774e-06,-5.97554e-10,4.4514e-14,-2795.77,-70.1137], Tmin=(890.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[CH]CC([O])C#C(22715)',
    structure = SMILES('[CH]CC([O])C#C'),
    E0 = (631.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,434.308,434.897,435.033,435.468,435.702],'cm^-1')),
        HinderedRotor(inertia=(0.000891727,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000892208,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3811,'amu*angstrom^2'), symmetry=1, barrier=(50.9815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979078,0.0564112,-5.52484e-05,2.72575e-08,-5.17409e-12,76059.3,24.8647], Tmin=(100,'K'), Tmax=(1400.37,'K')), NASAPolynomial(coeffs=[15.2218,0.0117614,-3.17255e-06,4.43032e-10,-2.59098e-14,72459.2,-47.2553], Tmin=(1400.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJ2_triplet) + radical(CC(C)OJ)"""),
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
    E0 = (219.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (304.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (341.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (341.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (301.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (416.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (219.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (449.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (373.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (480.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (361.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (466.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (328.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (333.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (410.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (331.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (308.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (497.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (641.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (554.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (714.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (451.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (669.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (615.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (768.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (433.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (462.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (277.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (350.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (244.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (244.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (627.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (513.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (483.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (767.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (787.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (227.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (420.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (531.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (715.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (941.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (295.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (395.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (219.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (401.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (395.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (420.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (452.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (376.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (463.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (298.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (308.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (441.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (227.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (699.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C1OC1CC=C[O](25778)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=[C]C1CC(C=O)O1(25421)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C1OC=CCC1[O](25509)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_O] for rate rule [R7_SMSS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#CC(=O)CC=C[O](25779)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59e+07,'cm^3/(mol*s)'), n=1.84, Ea=(32.6352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2830 used for CO-CtCs_O;HJ
Exact match found for rate rule [CO-CtCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C#CC([O])CC=C=O(25780)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.43214e-05,'m^3/(mol*s)'), n=3.00879, Ea=(45.2872,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-OneDeHH] for rate rule [CO-CtH_O;CsJ-CdHH]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond
Ea raised from 40.3 to 45.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=CCC=O(4687)', '[C]#C(5143)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CsH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C=C(O)CC=C[O](25781)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_Ct]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([O])CC=[C]O(25782)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(O)[CH]C=C[O](25783)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C[C]=CO(25784)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC(O)C[C]=C[O](25785)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC([O])[CH]C=CO(25786)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CC(O)CC=C[O](25787)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(O)C[CH][C]=O(25788)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5H_SSSD;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C=C([O])CC=CO(25789)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CC([O])CC=CO(25790)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C[O](602)', 'C#CC([CH2])[O](5294)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.00218e+08,'m^3/(mol*s)'), n=-0.446058, Ea=(0.74957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cs;Y_rad] for rate rule [C_rad/H2/Cs;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', 'C#CC([O])[CH]C=C[O](25791)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.66284e+07,'m^3/(mol*s)'), n=0.108445, Ea=(5.74998,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][CH]CC=C[O](1115)', '[C]#C(5143)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH]=C=C([O])CC=C[O](25792)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', 'C#CC([O])C[C]=C[O](25793)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', 'C#CC([O])C[CH][C]=O(25794)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[C]#CC([O])CC=C[O](25795)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC([O])CC1[CH]O1(25796)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[O]C=CCC1[C]=CO1(25704)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC1C[CH]C([O])O1(25797)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[O]C1[C]=COC=CC1(25563)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(=O)CC=CO(25798)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(O)CC=C=O(25799)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([O])[CH]C[CH][O](25800)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C=C([O])CC[CH][O](25801)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=C([O])C[CH]C[O](25802)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([O])CC[CH][O](25803)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]#CC([O])C[CH]C[O](25804)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC1CC=COO1(25617)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.214e+09,'s^-1'), n=0.749, Ea=(200.242,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for hex_1_ene_5_yne
Exact match found for rate rule [hex_1_ene_5_yne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C1[CH]CC([O])[C]=C1(25805)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O(T)(63)', 'C#C[CH]CC=C[O](20022)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(T)(63)', '[CH]=CCC([O])C#C(24341)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C1C([O])CC1C=O(25339)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(413769,'s^-1'), n=1.87624, Ea=(75.3751,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['H(8)', 'C#CC([O])C=CC=O(25806)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.76955,'m^3/(mol*s)'), n=1.94497, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDeH;HJ] for rate rule [Cds-CsH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(31.3381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 26.3 to 31.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#CC([O])[CH]CC=O(25486)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C#CC([O])CC[C]=O(22575)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[CH]=C=C([O])CCC=O(25807)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.07201e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[C]#CC([O])CCC=O(25808)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(13813.5,'s^-1'), n=1.88327, Ea=(38.7799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['[O]C1[C]=CC(C=O)C1(25389)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_csHCO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC1C[CH][CH]OO1(25601)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(243.719,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 241.3 to 243.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(O)C=CC=O(25809)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC(=O)CCC=O(25810)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([O])C([CH2])C=O(22594)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C#CC1CC(C=O)O1(22600)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=O(373)', '[CH]CC([O])C#C(22715)'],
    products = ['C#CC([O])CC=C[O](22586)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4887',
    isomers = [
        'C#CC([O])CC=C[O](22586)',
    ],
    reactants = [
        ('C=CC=O(5269)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4887',
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


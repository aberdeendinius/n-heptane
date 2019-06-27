species(
    label = 'C=CC([O])C=C=C[O](22584)',
    structure = SMILES('C=CC([O])C=C=C[O]'),
    E0 = (214.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11023,'amu*angstrom^2'), symmetry=1, barrier=(25.5263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10999,'amu*angstrom^2'), symmetry=1, barrier=(25.5209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438181,0.0654159,-3.65856e-05,-7.54109e-09,9.51967e-12,25941.3,30.6225], Tmin=(100,'K'), Tmax=(987.953,'K')), NASAPolynomial(coeffs=[19.5037,0.016784,-6.11062e-06,1.15527e-09,-8.53236e-14,20780.3,-68.1804], Tmin=(987.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C1OC1C=C=C[O](25839)',
    structure = SMILES('[CH2]C1OC1C=C=C[O]'),
    E0 = (248.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428816,0.0613023,-6.41526e-06,-5.80412e-08,3.46247e-11,30050.7,27.1792], Tmin=(100,'K'), Tmax=(879.081,'K')), NASAPolynomial(coeffs=[23.5523,0.00586335,3.24487e-06,-9.53752e-10,7.1335e-14,24061.9,-92.3376], Tmin=(879.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = 'C=CC1OC1[C]=C[O](25522)',
    structure = SMILES('C=CC1OC1[C]=C[O]'),
    E0 = (237.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788371,0.0502363,2.26384e-05,-8.60893e-08,4.38268e-11,28758.5,28.4947], Tmin=(100,'K'), Tmax=(890.6,'K')), NASAPolynomial(coeffs=[23.2305,0.0055508,3.3961e-06,-9.43326e-10,6.75579e-14,22535.8,-89.6685], Tmin=(890.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1OC=C=CC1[O](25422)',
    structure = SMILES('[CH2]C1OC=C=CC1[O]'),
    E0 = (255.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97877,0.128757,-0.000170309,1.07651e-07,-2.60956e-11,30920.2,14.8844], Tmin=(100,'K'), Tmax=(1022.22,'K')), NASAPolynomial(coeffs=[26.4157,0.017648,-7.26767e-06,1.31955e-09,-9.04972e-14,25115.1,-122.727], Tmin=(1022.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CJC(C)OC) + radical(CC(C)OJ)"""),
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
    label = 'C=CC(=O)C=C=C[O](25840)',
    structure = SMILES('C=CC(=O)C=C=C[O]'),
    E0 = (50.0932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05477,'amu*angstrom^2'), symmetry=1, barrier=(24.2512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05804,'amu*angstrom^2'), symmetry=1, barrier=(24.3264,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06919,0.0632461,-5.54715e-05,2.46332e-08,-4.44317e-12,6131.16,25.9385], Tmin=(100,'K'), Tmax=(1306.05,'K')), NASAPolynomial(coeffs=[13.5337,0.0250714,-1.16277e-05,2.2533e-09,-1.59268e-13,2875.32,-37.5236], Tmin=(1306.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.0932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([O])C=C=C=O(25841)',
    structure = SMILES('C=CC([O])C=C=C=O'),
    E0 = (258.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000736366,'amu*angstrom^2'), symmetry=1, barrier=(8.36075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71073,0.0418661,-1.79518e-05,-6.51091e-09,5.13273e-12,31194.8,13.2609], Tmin=(100,'K'), Tmax=(1072.44,'K')), NASAPolynomial(coeffs=[12.6875,0.0176199,-7.3901e-06,1.4394e-09,-1.04427e-13,27880.3,-44.9394], Tmin=(1072.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(64)',
    structure = SMILES('[CH]=C'),
    E0 = (289.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,826.012,826.012,3240.27],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90671,-0.00406241,3.8678e-05,-4.62976e-08,1.729e-11,34797.2,6.09789], Tmin=(100,'K'), Tmax=(931.962,'K')), NASAPolynomial(coeffs=[5.44797,0.00498356,-1.08821e-06,1.79837e-10,-1.45096e-14,33829.8,-4.87808], Tmin=(931.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C=CC=O(22476)',
    structure = SMILES('[O]C=C=CC=O'),
    E0 = (-8.09007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.29808,'amu*angstrom^2'), symmetry=1, barrier=(29.8455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63646,0.0447092,-3.72062e-05,1.25679e-08,-1.05205e-12,-881.798,19.0558], Tmin=(100,'K'), Tmax=(1155.06,'K')), NASAPolynomial(coeffs=[14.233,0.010485,-4.96566e-06,1.00356e-09,-7.36229e-14,-4418.66,-46.2447], Tmin=(1155.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.09007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=COJ)"""),
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
    label = '[CH2]C=C(O)C=C=C[O](25842)',
    structure = SMILES('[CH2]C=C(O)C=C=C[O]'),
    E0 = (35.3071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00381005,0.0707497,-3.00863e-05,-3.27647e-08,2.40608e-11,4406.37,27.0761], Tmin=(100,'K'), Tmax=(911.711,'K')), NASAPolynomial(coeffs=[25.2931,0.00574938,1.25219e-06,-3.96769e-10,2.59691e-14,-2114.76,-103.066], Tmin=(911.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.3071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC([O])C=C=[C]O(25843)',
    structure = SMILES('C=CC([O])C=C=[C]O'),
    E0 = (312.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.849186,'amu*angstrom^2'), symmetry=1, barrier=(19.5245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.850947,'amu*angstrom^2'), symmetry=1, barrier=(19.5649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852839,'amu*angstrom^2'), symmetry=1, barrier=(19.6084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2041,0.0745296,-6.82269e-05,2.67771e-08,-2.69425e-12,37766.6,33.1808], Tmin=(100,'K'), Tmax=(1017,'K')), NASAPolynomial(coeffs=[18.8073,0.0169008,-6.14927e-06,1.10872e-09,-7.78512e-14,33179,-60.8336], Tmin=(1017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]C(O)C=C=C[O](25844)',
    structure = SMILES('C=[C]C(O)C=C=C[O]'),
    E0 = (222.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.889253,'amu*angstrom^2'), symmetry=1, barrier=(20.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888425,'amu*angstrom^2'), symmetry=1, barrier=(20.4266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888295,'amu*angstrom^2'), symmetry=1, barrier=(20.4236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0160647,0.0765635,-7.56417e-05,3.67375e-08,-6.91216e-12,26854.1,32.81], Tmin=(100,'K'), Tmax=(1306.09,'K')), NASAPolynomial(coeffs=[19.9411,0.0155418,-5.5605e-06,9.66045e-10,-6.51257e-14,21649.3,-68.6373], Tmin=(1306.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=CC(O)[C]=C=C[O](25845)',
    structure = SMILES('C=CC(O)[C]=C=C[O]'),
    E0 = (222.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.889253,'amu*angstrom^2'), symmetry=1, barrier=(20.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888425,'amu*angstrom^2'), symmetry=1, barrier=(20.4266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888295,'amu*angstrom^2'), symmetry=1, barrier=(20.4236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0160647,0.0765635,-7.56417e-05,3.67375e-08,-6.91216e-12,26854.1,32.81], Tmin=(100,'K'), Tmax=(1306.09,'K')), NASAPolynomial(coeffs=[19.9411,0.0155418,-5.5605e-06,9.66045e-10,-6.51257e-14,21649.3,-68.6373], Tmin=(1306.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C=C=C[O](25846)',
    structure = SMILES('[CH]=CC(O)C=C=C[O]'),
    E0 = (231.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.947001,'amu*angstrom^2'), symmetry=1, barrier=(21.7734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.946066,'amu*angstrom^2'), symmetry=1, barrier=(21.7519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.947188,'amu*angstrom^2'), symmetry=1, barrier=(21.7777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.196495,0.0729953,-5.93798e-05,1.51892e-08,1.85962e-12,27960.7,32.0848], Tmin=(100,'K'), Tmax=(995.154,'K')), NASAPolynomial(coeffs=[19.8353,0.0156986,-5.63582e-06,1.03763e-09,-7.49098e-14,22980.4,-67.9496], Tmin=(995.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([O])C#C[CH]O(25847)',
    structure = SMILES('C=CC([O])C#C[CH]O'),
    E0 = (292.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2100,2250,500,550,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,266.49,266.49,266.499,266.503],'cm^-1')),
        HinderedRotor(inertia=(0.00237366,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00237376,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65688,'amu*angstrom^2'), symmetry=1, barrier=(83.5014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65692,'amu*angstrom^2'), symmetry=1, barrier=(83.5011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481179,0.0749326,-8.33407e-05,4.98302e-08,-1.18066e-11,35273.4,32.6797], Tmin=(100,'K'), Tmax=(1034.24,'K')), NASAPolynomial(coeffs=[13.8893,0.0230746,-8.12732e-06,1.34686e-09,-8.6749e-14,32500,-32.4578], Tmin=(1034.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CC(C)OJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=CC(O)[CH][C]=C=O(25848)',
    structure = SMILES('C=CC(O)[CH][C]=C=O'),
    E0 = (78.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707426,0.0578897,-1.66437e-05,-2.60992e-08,1.54472e-11,9531.47,31.3698], Tmin=(100,'K'), Tmax=(988.137,'K')), NASAPolynomial(coeffs=[18.715,0.0182173,-6.85302e-06,1.32057e-09,-9.85008e-14,4350.72,-63.4988], Tmin=(988.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCCJ=C=O) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]C=C([O])C=C=CO(25849)',
    structure = SMILES('[CH2]C=C([O])C=C=CO'),
    E0 = (31.6493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37449,0.0892768,-9.46512e-05,4.71129e-08,-8.56681e-12,4024.96,32.127], Tmin=(100,'K'), Tmax=(1610.77,'K')), NASAPolynomial(coeffs=[24.2289,0.00557578,2.0307e-06,-6.5675e-10,5.09003e-14,-1613.09,-95.4977], Tmin=(1610.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.6493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C([O])C=C=CO(25850)',
    structure = SMILES('C=[C]C([O])C=C=CO'),
    E0 = (310.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.992004,'amu*angstrom^2'), symmetry=1, barrier=(22.8081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993424,'amu*angstrom^2'), symmetry=1, barrier=(22.8408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992967,'amu*angstrom^2'), symmetry=1, barrier=(22.8303,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00116786,0.0769443,-6.47385e-05,1.56492e-08,3.1423e-12,37547.3,31.262], Tmin=(100,'K'), Tmax=(961.26,'K')), NASAPolynomial(coeffs=[21.4783,0.0127431,-3.84676e-06,6.6847e-10,-4.8602e-14,32254.5,-77.5665], Tmin=(961.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C=C=CO(25851)',
    structure = SMILES('[CH]=CC([O])C=C=CO'),
    E0 = (320.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05943,'amu*angstrom^2'), symmetry=1, barrier=(24.3584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05787,'amu*angstrom^2'), symmetry=1, barrier=(24.3226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06038,'amu*angstrom^2'), symmetry=1, barrier=(24.3802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.848316,0.0851095,-8.77431e-05,4.25486e-08,-7.70223e-12,38699.3,34.2493], Tmin=(100,'K'), Tmax=(1537.27,'K')), NASAPolynomial(coeffs=[24.6427,0.00737192,-7.57372e-07,-2.40118e-12,2.82857e-15,32210.1,-95.3069], Tmin=(1537.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=[C]C=C[O](23191)',
    structure = SMILES('[O]C=[C]C=C[O]'),
    E0 = (158.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62308,'amu*angstrom^2'), symmetry=1, barrier=(37.3178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34928,0.0410469,6.39572e-06,-6.06125e-08,3.33202e-11,19145.6,19.648], Tmin=(100,'K'), Tmax=(893.363,'K')), NASAPolynomial(coeffs=[22.6723,-0.00687368,7.01377e-06,-1.49144e-09,1.02026e-13,13438.2,-91.4392], Tmin=(893.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[CH2]C=C([O])C=C=C[O](25852)',
    structure = SMILES('[CH2]C=C([O])C=C=C[O]'),
    E0 = (173.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.56956,'amu*angstrom^2'), symmetry=1, barrier=(36.0873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56317,'amu*angstrom^2'), symmetry=1, barrier=(35.9404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.708747,0.079265,-8.03759e-05,3.88368e-08,-6.93886e-12,21011,31.1506], Tmin=(100,'K'), Tmax=(1612.54,'K')), NASAPolynomial(coeffs=[21.9305,0.00796422,3.49031e-08,-2.30799e-10,2.08404e-14,15678.5,-82.7833], Tmin=(1612.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C([O])C=C=C[O](25853)',
    structure = SMILES('C=[C]C([O])C=C=C[O]'),
    E0 = (452.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07032,'amu*angstrom^2'), symmetry=1, barrier=(24.6088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07258,'amu*angstrom^2'), symmetry=1, barrier=(24.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374667,0.0703608,-6.24987e-05,2.30454e-08,-1.93579e-12,54545.9,31.3249], Tmin=(100,'K'), Tmax=(1041.4,'K')), NASAPolynomial(coeffs=[18.5633,0.0158704,-6.15375e-06,1.14956e-09,-8.21564e-14,49924.1,-61.1646], Tmin=(1041.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=CC([O])[C]=C=C[O](25854)',
    structure = SMILES('C=CC([O])[C]=C=C[O]'),
    E0 = (452.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07032,'amu*angstrom^2'), symmetry=1, barrier=(24.6088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07258,'amu*angstrom^2'), symmetry=1, barrier=(24.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374667,0.0703608,-6.24987e-05,2.30454e-08,-1.93579e-12,54545.9,31.3249], Tmin=(100,'K'), Tmax=(1041.4,'K')), NASAPolynomial(coeffs=[18.5633,0.0158704,-6.15375e-06,1.14956e-09,-8.21564e-14,49924.1,-61.1646], Tmin=(1041.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C=C=C[O](25855)',
    structure = SMILES('[CH]=CC([O])C=C=C[O]'),
    E0 = (461.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13381,'amu*angstrom^2'), symmetry=1, barrier=(26.0686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13401,'amu*angstrom^2'), symmetry=1, barrier=(26.0732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360662,0.068984,-5.34068e-05,1.00603e-08,3.5076e-12,55661.2,31.3038], Tmin=(100,'K'), Tmax=(993.696,'K')), NASAPolynomial(coeffs=[19.7627,0.0139192,-5.05812e-06,9.52127e-10,-7.0093e-14,50667.9,-67.8996], Tmin=(993.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC([O])[CH][C]=C=O(25856)',
    structure = SMILES('C=CC([O])[CH][C]=C=O'),
    E0 = (308.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,463.136,464.672,465.432,466.061],'cm^-1')),
        HinderedRotor(inertia=(0.226813,'amu*angstrom^2'), symmetry=1, barrier=(34.6169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226299,'amu*angstrom^2'), symmetry=1, barrier=(34.6341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224829,'amu*angstrom^2'), symmetry=1, barrier=(34.6039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870503,0.0538915,-1.07165e-05,-3.11685e-08,1.70699e-11,37232,30.5927], Tmin=(100,'K'), Tmax=(987.65,'K')), NASAPolynomial(coeffs=[18.6473,0.0164297,-6.27061e-06,1.23396e-09,-9.35921e-14,32036.2,-63.476], Tmin=(987.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CCCJ=C=O) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC([O])C=C1[CH]O1(25857)',
    structure = SMILES('C=CC([O])C=C1[CH]O1'),
    E0 = (240.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850301,0.0443357,4.02116e-05,-9.93639e-08,4.51491e-11,29086.8,27.7241], Tmin=(100,'K'), Tmax=(942.625,'K')), NASAPolynomial(coeffs=[24.5617,0.00669972,-1.22758e-07,4.56836e-11,-1.54401e-14,21818.4,-100.111], Tmin=(942.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CC(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = '[O]C=C=CC1[CH]CO1(25858)',
    structure = SMILES('[O]C=C=CC1[CH]CO1'),
    E0 = (244.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17943,0.0401437,4.73873e-05,-1.08238e-07,5.08523e-11,29477.5,27.7526], Tmin=(100,'K'), Tmax=(893.452,'K')), NASAPolynomial(coeffs=[21.6737,0.00755493,2.77022e-06,-8.29399e-10,5.90257e-14,23454,-82.0264], Tmin=(893.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Oxetane) + radical(CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[CH2][CH]C1C=C(C=O)O1(25508)',
    structure = SMILES('[CH2][CH]C1C=C(C=O)O1'),
    E0 = (222.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707127,0.0521251,9.2031e-06,-5.92353e-08,2.8293e-11,26911.6,28.0093], Tmin=(100,'K'), Tmax=(978.653,'K')), NASAPolynomial(coeffs=[22.3464,0.0123878,-4.54608e-06,9.86521e-10,-8.20875e-14,20343.6,-87.8376], Tmin=(978.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC([O])C1[C]=CO1(25859)',
    structure = SMILES('C=CC([O])C1[C]=CO1'),
    E0 = (338.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632794,0.047374,3.95339e-05,-1.04428e-07,4.83955e-11,40812.3,28.3745], Tmin=(100,'K'), Tmax=(937.078,'K')), NASAPolynomial(coeffs=[26.9029,0.00308634,1.81959e-06,-3.30788e-10,1.02228e-14,32910,-112.551], Tmin=(937.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC1C=[C]C([O])O1(25652)',
    structure = SMILES('C=CC1C=[C]C([O])O1'),
    E0 = (234.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31949,0.0497032,-1.78742e-05,-6.73069e-09,4.23751e-12,28364.1,26.8646], Tmin=(100,'K'), Tmax=(1200.65,'K')), NASAPolynomial(coeffs=[12.8479,0.0285402,-1.29783e-05,2.51302e-09,-1.77995e-13,24352.9,-36.0376], Tmin=(1200.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(25dihydrofuran) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[CH]COC=C=C1(25753)',
    structure = SMILES('[O]C1[CH]COC=C=C1'),
    E0 = (308.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03508,0.0416145,4.3683e-05,-9.78678e-08,4.32773e-11,37252.5,21.6092], Tmin=(100,'K'), Tmax=(947.987,'K')), NASAPolynomial(coeffs=[22.2191,0.0114827,-2.39627e-06,4.71109e-10,-4.41581e-14,30573.6,-93.5028], Tmin=(947.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(=O)C=C=CO(25860)',
    structure = SMILES('C=CC(=O)C=C=CO'),
    E0 = (-91.3694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355446,0.0737493,-7.11446e-05,3.43182e-08,-6.52396e-12,-10852.6,27.0925], Tmin=(100,'K'), Tmax=(1276.29,'K')), NASAPolynomial(coeffs=[17.3765,0.0204039,-8.44859e-06,1.56905e-09,-1.09045e-13,-15197.4,-59.1765], Tmin=(1276.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.3694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(O)C=C=C=O(25861)',
    structure = SMILES('C=CC(O)C=C=C=O'),
    E0 = (28.2618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53767,0.0459755,-2.42338e-05,-1.03653e-09,3.36334e-12,3494.7,14.0742], Tmin=(100,'K'), Tmax=(1087.76,'K')), NASAPolynomial(coeffs=[12.8423,0.0192655,-7.89304e-06,1.50765e-09,-1.07838e-13,156.201,-45.4562], Tmin=(1087.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.2618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]CC([O])=C[C]=C[O](25862)',
    structure = SMILES('[CH2]CC([O])=C[C]=C[O]'),
    E0 = (295.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,252.493,252.524,252.535,252.558],'cm^-1')),
        HinderedRotor(inertia=(0.470626,'amu*angstrom^2'), symmetry=1, barrier=(21.2928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470739,'amu*angstrom^2'), symmetry=1, barrier=(21.2953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470657,'amu*angstrom^2'), symmetry=1, barrier=(21.2942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.741056,0.0849086,-9.01881e-05,4.61505e-08,-8.80114e-12,35716.7,33.944], Tmin=(100,'K'), Tmax=(1479.78,'K')), NASAPolynomial(coeffs=[22.4853,0.00924317,-4.3053e-07,-1.69873e-10,1.82997e-14,30253.2,-82.446], Tmin=(1479.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJ) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=C([O])C[C]=C[O](25863)',
    structure = SMILES('[CH2]C=C([O])C[C]=C[O]'),
    E0 = (297.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,385.468,385.499,385.545,385.685],'cm^-1')),
        HinderedRotor(inertia=(0.171114,'amu*angstrom^2'), symmetry=1, barrier=(18.0275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170783,'amu*angstrom^2'), symmetry=1, barrier=(18.0288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170892,'amu*angstrom^2'), symmetry=1, barrier=(18.0258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286024,0.0694228,-4.46492e-05,-4.63873e-09,1.051e-11,35983.8,31.5419], Tmin=(100,'K'), Tmax=(939.532,'K')), NASAPolynomial(coeffs=[20.6509,0.0131447,-3.37266e-06,5.39049e-10,-3.89268e-14,30814.3,-72.583], Tmin=(939.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC([O])[C]C=C[O](25864)',
    structure = SMILES('C=CC([O])[C]C=C[O]'),
    E0 = (502.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300.752,300.752,300.752,300.752,300.752,300.753],'cm^-1')),
        HinderedRotor(inertia=(0.369974,'amu*angstrom^2'), symmetry=1, barrier=(23.7474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369973,'amu*angstrom^2'), symmetry=1, barrier=(23.7474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369973,'amu*angstrom^2'), symmetry=1, barrier=(23.7474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.306001,0.0662585,-3.13156e-05,-1.88695e-08,1.50764e-11,60540,32.5495], Tmin=(100,'K'), Tmax=(961.058,'K')), NASAPolynomial(coeffs=[21.6348,0.0131125,-3.9709e-06,7.33796e-10,-5.66335e-14,54795,-78.0622], Tmin=(961.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC([O])[CH]C=[C][O](25865)',
    structure = SMILES('C=CC([O])[CH]C=[C][O]'),
    E0 = (363.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,529.299,529.3,529.346,529.359,529.444],'cm^-1')),
        HinderedRotor(inertia=(0.138617,'amu*angstrom^2'), symmetry=1, barrier=(27.5615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138611,'amu*angstrom^2'), symmetry=1, barrier=(27.5622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138604,'amu*angstrom^2'), symmetry=1, barrier=(27.5626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826801,0.0507406,1.14793e-05,-6.03735e-08,2.88512e-11,43887.9,34.7022], Tmin=(100,'K'), Tmax=(962.982,'K')), NASAPolynomial(coeffs=[20.7236,0.0145845,-4.61818e-06,9.04133e-10,-7.22912e-14,37900.3,-71.7302], Tmin=(962.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(C=CJO) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])C[C]=C[O](25866)',
    structure = SMILES('C=[C]C([O])C[C]=C[O]'),
    E0 = (529.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,385.451,385.886,386.378,386.661,386.978],'cm^-1')),
        HinderedRotor(inertia=(0.159193,'amu*angstrom^2'), symmetry=1, barrier=(16.7341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159137,'amu*angstrom^2'), symmetry=1, barrier=(16.7346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159172,'amu*angstrom^2'), symmetry=1, barrier=(16.7276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.254372,0.0728415,-6.28551e-05,2.18098e-08,-1.20139e-12,63840.3,33.7772], Tmin=(100,'K'), Tmax=(1023.62,'K')), NASAPolynomial(coeffs=[18.6718,0.0176677,-6.61653e-06,1.2123e-09,-8.58081e-14,59189.9,-59.8043], Tmin=(1023.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2]CC([O])[C]=C=C[O](25867)',
    structure = SMILES('[CH2]CC([O])[C]=C=C[O]'),
    E0 = (533.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,387.727,389.064,389.633,389.64],'cm^-1')),
        HinderedRotor(inertia=(0.167926,'amu*angstrom^2'), symmetry=1, barrier=(18.0512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169333,'amu*angstrom^2'), symmetry=1, barrier=(18.0563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16878,'amu*angstrom^2'), symmetry=1, barrier=(18.0442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0932282,0.0788399,-7.90074e-05,3.88757e-08,-7.401e-12,64332.7,34.4847], Tmin=(100,'K'), Tmax=(1292.39,'K')), NASAPolynomial(coeffs=[20.3764,0.015485,-5.47441e-06,9.44115e-10,-6.34394e-14,59041.9,-69.5195], Tmin=(1292.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(RCCJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=C([O])[CH]C=C[O](25868)',
    structure = SMILES('[CH2]C=C([O])[CH]C=C[O]'),
    E0 = (177.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,186.661,186.741,187.185,188.046],'cm^-1')),
        HinderedRotor(inertia=(1.34399,'amu*angstrom^2'), symmetry=1, barrier=(33.3606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34452,'amu*angstrom^2'), symmetry=1, barrier=(33.3193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32889,'amu*angstrom^2'), symmetry=1, barrier=(33.3202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.100206,0.0687762,-2.75559e-05,-3.02297e-08,2.12057e-11,21451,28.9347], Tmin=(100,'K'), Tmax=(935.727,'K')), NASAPolynomial(coeffs=[23.9777,0.00961323,-1.49749e-06,2.08805e-10,-1.91129e-14,15104,-94.7116], Tmin=(935.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C([O])[CH]C=C[O](25869)',
    structure = SMILES('C=[C]C([O])[CH]C=C[O]'),
    E0 = (361.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,400.605,401.786,401.917,403.537,404.217],'cm^-1')),
        HinderedRotor(inertia=(0.265475,'amu*angstrom^2'), symmetry=1, barrier=(30.3258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267804,'amu*angstrom^2'), symmetry=1, barrier=(30.3472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261097,'amu*angstrom^2'), symmetry=1, barrier=(30.3328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588368,0.053539,1.36759e-05,-6.99299e-08,3.40763e-11,43670.1,32.9026], Tmin=(100,'K'), Tmax=(949.779,'K')), NASAPolynomial(coeffs=[23.6005,0.0100815,-2.11859e-06,4.17679e-10,-3.92308e-14,36887.6,-89.6254], Tmin=(949.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]C([O])[C]=C=C[O](25870)',
    structure = SMILES('C[CH]C([O])[C]=C=C[O]'),
    E0 = (528.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1685,370,3010,987.5,1337.5,450,1655,418.679,418.684,418.685,418.689],'cm^-1')),
        HinderedRotor(inertia=(0.138605,'amu*angstrom^2'), symmetry=1, barrier=(17.2424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13861,'amu*angstrom^2'), symmetry=1, barrier=(17.2425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138629,'amu*angstrom^2'), symmetry=1, barrier=(17.2425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291736,0.0703216,-5.15119e-05,7.10909e-09,4.74567e-12,63676.7,33.8063], Tmin=(100,'K'), Tmax=(983.335,'K')), NASAPolynomial(coeffs=[19.5974,0.0160058,-5.59644e-06,1.02343e-09,-7.4078e-14,58709.1,-64.9611], Tmin=(983.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CCJCO) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C([O])C=[C]C[O](25871)',
    structure = SMILES('[CH2]C=C([O])C=[C]C[O]'),
    E0 = (374.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,275.689,277.156,279.39,1771],'cm^-1')),
        HinderedRotor(inertia=(0.0807691,'amu*angstrom^2'), symmetry=1, barrier=(4.34892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813478,'amu*angstrom^2'), symmetry=1, barrier=(4.35244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2154,'amu*angstrom^2'), symmetry=1, barrier=(64.5975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10619,0.0696283,-7.29382e-05,3.56711e-08,-1.67674e-12,45110.4,29.4312], Tmin=(100,'K'), Tmax=(623.393,'K')), NASAPolynomial(coeffs=[8.40916,0.0341229,-1.48256e-05,2.74119e-09,-1.87582e-13,43979.3,-4.11972], Tmin=(623.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(CCOJ) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C[C]=C[O](25872)',
    structure = SMILES('[CH]=CC([O])C[C]=C[O]'),
    E0 = (538.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,326.789,326.809,326.838,326.891],'cm^-1')),
        HinderedRotor(inertia=(0.239658,'amu*angstrom^2'), symmetry=1, barrier=(18.1676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239694,'amu*angstrom^2'), symmetry=1, barrier=(18.1679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239666,'amu*angstrom^2'), symmetry=1, barrier=(18.1682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.236675,0.0715017,-5.38578e-05,8.89105e-09,4.24093e-12,64955.7,33.7697], Tmin=(100,'K'), Tmax=(984.54,'K')), NASAPolynomial(coeffs=[19.923,0.0156322,-5.47389e-06,1.00403e-09,-7.28626e-14,59910.7,-66.8334], Tmin=(984.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])C=[C]C[O](25873)',
    structure = SMILES('C=[C]C([O])C=[C]C[O]'),
    E0 = (653.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,295.483,295.49,295.5,1747.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.13353,'amu*angstrom^2'), symmetry=1, barrier=(8.27397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133534,'amu*angstrom^2'), symmetry=1, barrier=(8.27397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52524,'amu*angstrom^2'), symmetry=1, barrier=(32.5446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.915659,0.0763685,-0.000113977,1.05573e-07,-3.94399e-11,78701.2,34.1457], Tmin=(100,'K'), Tmax=(802.227,'K')), NASAPolynomial(coeffs=[4.67707,0.0420859,-2.08422e-05,4.04872e-09,-2.82671e-13,78597.4,19.9422], Tmin=(802.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([O])[CH]C=C[O](25874)',
    structure = SMILES('[CH]=CC([O])[CH]C=C[O]'),
    E0 = (371.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,341.855,341.864,341.879,341.904],'cm^-1')),
        HinderedRotor(inertia=(0.37829,'amu*angstrom^2'), symmetry=1, barrier=(31.3713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378218,'amu*angstrom^2'), symmetry=1, barrier=(31.3711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378314,'amu*angstrom^2'), symmetry=1, barrier=(31.3718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.545532,0.0524831,2.17545e-05,-8.17954e-08,3.91434e-11,44786.6,32.986], Tmin=(100,'K'), Tmax=(944.772,'K')), NASAPolynomial(coeffs=[25.0404,0.00773179,-7.97378e-07,1.67682e-10,-2.28529e-14,37527,-97.7208], Tmin=(944.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ) + radical(C=COJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]CC([O])[CH][C]=C=O(25875)',
    structure = SMILES('[CH2]CC([O])[CH][C]=C=O'),
    E0 = (433.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,180,180,270.528,737.103],'cm^-1')),
        HinderedRotor(inertia=(0.0586692,'amu*angstrom^2'), symmetry=1, barrier=(22.6202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269145,'amu*angstrom^2'), symmetry=1, barrier=(103.765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586753,'amu*angstrom^2'), symmetry=1, barrier=(22.6204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.75398,'amu*angstrom^2'), symmetry=1, barrier=(63.3195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483137,0.0788692,-8.32803e-05,4.63177e-08,-1.04458e-11,52214.2,28.7946], Tmin=(100,'K'), Tmax=(1065.5,'K')), NASAPolynomial(coeffs=[13.7131,0.029203,-1.33612e-05,2.57076e-09,-1.81472e-13,49394.9,-35.8716], Tmin=(1065.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJCO) + radical(RCCJ) + radical(CC(C)OJ) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH]=CC([O])C=[C]C[O](25876)',
    structure = SMILES('[CH]=CC([O])C=[C]C[O]'),
    E0 = (662.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,281.583,282.389,1242.19,1972.94],'cm^-1')),
        HinderedRotor(inertia=(0.175404,'amu*angstrom^2'), symmetry=1, barrier=(9.92944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176935,'amu*angstrom^2'), symmetry=1, barrier=(9.92832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176726,'amu*angstrom^2'), symmetry=1, barrier=(9.93389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907094,0.0748788,-0.000104258,9.14539e-08,-3.33908e-11,79816.3,34.1082], Tmin=(100,'K'), Tmax=(769.497,'K')), NASAPolynomial(coeffs=[5.99245,0.0399613,-1.96565e-05,3.83176e-09,-2.69098e-13,79284.8,12.5386], Tmin=(769.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]C([O])[CH][C]=C=O(25877)',
    structure = SMILES('C[CH]C([O])[CH][C]=C=O'),
    E0 = (427.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2120,512.5,787.5,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,426.537,426.572,426.575,426.58],'cm^-1')),
        HinderedRotor(inertia=(0.000926804,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272387,'amu*angstrom^2'), symmetry=1, barrier=(35.1671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272344,'amu*angstrom^2'), symmetry=1, barrier=(35.1668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272252,'amu*angstrom^2'), symmetry=1, barrier=(35.1668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441027,0.0753468,-7.2957e-05,3.61944e-08,-7.19319e-12,51576.5,29.6496], Tmin=(100,'K'), Tmax=(1210.11,'K')), NASAPolynomial(coeffs=[15.5542,0.0253905,-1.10331e-05,2.07957e-09,-1.45299e-13,47918.8,-46.1447], Tmin=(1210.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCCJ=C=O) + radical(CCJCO) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4427.86,'J/mol'), sigma=(7.0864,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.62 K, Pc=28.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349841,0.0658846,-3.00643e-05,-2.26308e-08,1.7727e-11,26589.5,30.1039], Tmin=(100,'K'), Tmax=(927.933,'K')), NASAPolynomial(coeffs=[21.5129,0.0116462,-2.17941e-06,2.92247e-10,-2.20858e-14,21069.4,-78.9933], Tmin=(927.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=C1[CH]C([O])[CH]C1(25878)',
    structure = SMILES('[O]C=C1[CH]C([O])[CH]C1'),
    E0 = (300.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11698,0.0353796,7.07707e-05,-1.29698e-07,5.52128e-11,36273,25.155], Tmin=(100,'K'), Tmax=(946.172,'K')), NASAPolynomial(coeffs=[23.8122,0.0102211,-1.56506e-06,3.3862e-10,-3.79873e-14,28809.7,-99.8245], Tmin=(946.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclopentane) + radical(CCJCO) + radical(C=COJ) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
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
    label = '[CH2]C=CC=C=C[O](23648)',
    structure = SMILES('[CH2]C=CC=C=C[O]'),
    E0 = (249.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.60533,'amu*angstrom^2'), symmetry=1, barrier=(36.9096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59987,'amu*angstrom^2'), symmetry=1, barrier=(36.7841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0576,0.0479059,8.36085e-06,-5.70423e-08,2.88683e-11,30093.3,24.0495], Tmin=(100,'K'), Tmax=(928.436,'K')), NASAPolynomial(coeffs=[19.6188,0.0110491,-1.74335e-06,2.26068e-10,-1.93777e-14,24788.7,-74.1257], Tmin=(928.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=CC([O])C=C(20057)',
    structure = SMILES('[CH]=C=CC([O])C=C'),
    E0 = (436.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.970908,'amu*angstrom^2'), symmetry=1, barrier=(22.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9706,'amu*angstrom^2'), symmetry=1, barrier=(22.316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0173,0.0582997,-4.56264e-05,1.57995e-08,-1.35038e-12,52592.3,27.902], Tmin=(100,'K'), Tmax=(1075.48,'K')), NASAPolynomial(coeffs=[13.8762,0.0208721,-7.92724e-06,1.42018e-09,-9.74795e-14,49225.1,-37.8661], Tmin=(1075.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C(C=O)=CC1[O](25338)',
    structure = SMILES('[CH2]C1C(C=O)=CC1[O]'),
    E0 = (244.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09601,0.0500875,-4.52294e-06,-3.13859e-08,1.54828e-11,29528.7,27.8205], Tmin=(100,'K'), Tmax=(1014.32,'K')), NASAPolynomial(coeffs=[16.0862,0.0220447,-9.00186e-06,1.75797e-09,-1.29547e-13,24889.3,-52.591], Tmin=(1014.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC([O])C#CC=O(25879)',
    structure = SMILES('C=CC([O])C#CC=O'),
    E0 = (177.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2100,2250,500,550,3010,987.5,1337.5,450,1655,356.319,356.324,356.342],'cm^-1')),
        HinderedRotor(inertia=(0.228109,'amu*angstrom^2'), symmetry=1, barrier=(20.5509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22808,'amu*angstrom^2'), symmetry=1, barrier=(20.5508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401202,'amu*angstrom^2'), symmetry=1, barrier=(36.1498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975915,0.0617986,-5.36149e-05,2.38828e-08,-4.28146e-12,21433,29.0097], Tmin=(100,'K'), Tmax=(1330.2,'K')), NASAPolynomial(coeffs=[14.0071,0.0226137,-9.42886e-06,1.73806e-09,-1.19613e-13,17966.1,-37.5764], Tmin=(1330.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC([O])[C]=CC=O(25880)',
    structure = SMILES('C=CC([O])[C]=CC=O'),
    E0 = (255.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,331.144,331.149,331.167],'cm^-1')),
        HinderedRotor(inertia=(0.169585,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169585,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169587,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42684,0.0636407,-5.50752e-05,2.53375e-08,-5.05754e-12,30832,28.5063], Tmin=(100,'K'), Tmax=(1127.56,'K')), NASAPolynomial(coeffs=[9.1639,0.0361935,-1.85618e-05,3.74895e-09,-2.70949e-13,29087.2,-9.74941], Tmin=(1127.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])[CH]C=C=O(25881)',
    structure = SMILES('C=CC([O])[CH]C=C=O'),
    E0 = (106.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.887913,0.0528099,-3.28189e-06,-3.83465e-08,1.92745e-11,12903.2,30.8213], Tmin=(100,'K'), Tmax=(990.799,'K')), NASAPolynomial(coeffs=[18.2209,0.0191812,-7.39757e-06,1.44824e-09,-1.0885e-13,7684.41,-61.6434], Tmin=(990.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=C([O])C=CC=O(25882)',
    structure = SMILES('[CH2]C=C([O])C=CC=O'),
    E0 = (-23.6248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672595,0.0690087,-6.19815e-05,2.85352e-08,-5.28123e-12,-2718.21,27.1249], Tmin=(100,'K'), Tmax=(1291.18,'K')), NASAPolynomial(coeffs=[14.9713,0.0247113,-1.05188e-05,1.96316e-09,-1.36211e-13,-6410.56,-45.5117], Tmin=(1291.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.6248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])C=CC=O(25883)',
    structure = SMILES('C=[C]C([O])C=CC=O'),
    E0 = (255.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,331.149,331.15,331.151],'cm^-1')),
        HinderedRotor(inertia=(0.169584,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169581,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169581,'amu*angstrom^2'), symmetry=1, barrier=(13.1965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42684,0.0636407,-5.50752e-05,2.53375e-08,-5.05754e-12,30832,28.5063], Tmin=(100,'K'), Tmax=(1127.56,'K')), NASAPolynomial(coeffs=[9.1639,0.0361935,-1.85618e-05,3.74895e-09,-2.70949e-13,29087.2,-9.74941], Tmin=(1127.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C=CC=O(25884)',
    structure = SMILES('[CH]=CC([O])C=CC=O'),
    E0 = (264.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,241.063,245.621],'cm^-1')),
        HinderedRotor(inertia=(0.402029,'amu*angstrom^2'), symmetry=1, barrier=(14.9651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345689,'amu*angstrom^2'), symmetry=1, barrier=(14.887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337269,'amu*angstrom^2'), symmetry=1, barrier=(14.931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26288,0.0641803,-5.31865e-05,2.20801e-08,-3.81796e-12,31953.2,29.0106], Tmin=(100,'K'), Tmax=(1316.78,'K')), NASAPolynomial(coeffs=[12.0587,0.0313856,-1.58285e-05,3.16621e-09,-2.27018e-13,29110,-26.0439], Tmin=(1316.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[CH]CC(C=O)=C1(25670)',
    structure = SMILES('[O]C1[CH]CC(C=O)=C1'),
    E0 = (148.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52074,0.0371962,3.03895e-05,-6.45157e-08,2.59944e-11,17982.6,27.0711], Tmin=(100,'K'), Tmax=(1019.64,'K')), NASAPolynomial(coeffs=[15.7066,0.023183,-1.02492e-05,2.10457e-09,-1.59509e-13,12925.2,-52.2573], Tmin=(1019.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentene) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC1[CH][C]=COO1(25602)',
    structure = SMILES('C=CC1[CH][C]=COO1'),
    E0 = (337.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46508,0.0292424,7.92369e-05,-1.35483e-07,5.73933e-11,40701.2,24.2453], Tmin=(100,'K'), Tmax=(931.093,'K')), NASAPolynomial(coeffs=[21.4793,0.0110962,-8.12497e-07,7.99e-11,-1.50659e-14,34033.8,-86.6729], Tmin=(931.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(34dihydro12dioxin) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(=O)C=CC=O(25885)',
    structure = SMILES('C=CC(=O)C=CC=O'),
    E0 = (-146.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94543,0.0486496,-2.69518e-05,5.87281e-09,-4.6868e-13,-17626.5,20.0072], Tmin=(100,'K'), Tmax=(3031.82,'K')), NASAPolynomial(coeffs=[37.5942,0.0029364,-4.33526e-06,8.99687e-10,-5.86055e-14,-38636.5,-185.585], Tmin=(3031.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CC1C=C(C=O)O1(25505)',
    structure = SMILES('C=CC1C=C(C=O)O1'),
    E0 = (-58.5043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.781166,0.0511062,9.74568e-06,-5.68446e-08,2.6545e-11,-6903.02,23.7729], Tmin=(100,'K'), Tmax=(991.218,'K')), NASAPolynomial(coeffs=[21.335,0.0148517,-6.04474e-06,1.29562e-09,-1.04092e-13,-13271.3,-86.7758], Tmin=(991.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.5043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    label = '[C]=CC([O])C=C(13789)',
    structure = SMILES('[C]=CC([O])C=C'),
    E0 = (699.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,433.301,433.306,433.319],'cm^-1')),
        HinderedRotor(inertia=(0.101337,'amu*angstrom^2'), symmetry=1, barrier=(13.5023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101341,'amu*angstrom^2'), symmetry=1, barrier=(13.5024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78515,0.0459501,-3.62672e-05,1.43661e-08,-2.30688e-12,84197.7,24.0389], Tmin=(100,'K'), Tmax=(1454.69,'K')), NASAPolynomial(coeffs=[11.6236,0.0188971,-8.37159e-06,1.58187e-09,-1.09821e-13,81335.3,-27.1133], Tmin=(1454.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(CC(C)OJ)"""),
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
    E0 = (214.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (322.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (298.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (402.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (283.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (473.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (339.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (214.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (339.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (490.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (363.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (363.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (275.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (607.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (305.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (305.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (462.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (471.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (447.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (384.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (664.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (664.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (673.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (520.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (401.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (347.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (340.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (340.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (284.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (308.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (239.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (239.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (318.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (325.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (524.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (386.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (622.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (622.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (266.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (401.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (567.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (413.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (556.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (692.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (396.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (458.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (687.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (452.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (420.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (450.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (656.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (843.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (268.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (398.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (214.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (452.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (377.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (391.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (649.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (297.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (252.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (337.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (303.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (222.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (767.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C1OC1C=C=C[O](25839)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC1OC1[C]=C[O](25522)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C1OC=C=CC1[O](25422)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7_SMMS;doublebond_intra_2H_pri;radadd_intra] for rate rule [R7_SMMS;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C=CC(=O)C=C=C[O](25840)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C=CC([O])C=C=C=O(25841)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(64)', '[O]C=C=CC=O(22476)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdH_O;CdsJ-H]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(25.9901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;CdsJ=Cdd]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 20.7 to 26.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C=C(O)C=C=C[O](25842)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=CC([O])C=C=[C]O(25843)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=CC(O)[C]=C=C[O](25845)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(O)C=C=C[O](25846)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=CC([O])C#C[CH]O(25847)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC(O)[CH][C]=C=O(25848)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C=C([O])C=C=CO(25849)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C([O])C=C=CO(25850)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R6H;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC([O])C=C=CO(25851)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(64)', '[O]C=[C]C=C[O](23191)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C=C([O])C=C=C[O](25852)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C=[C]C([O])C=C=C[O](25853)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', 'C=CC([O])[C]=C=C[O](25854)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH]=CC([O])C=C=C[O](25855)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', 'C=CC([O])[CH][C]=C=O(25856)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC([O])C=C1[CH]O1(25857)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[O]C=C=CC1[CH]CO1(25858)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.15968e+08,'s^-1'), n=1.10215, Ea=(132.51,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2][CH]C1C=C(C=O)O1(25508)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC([O])C1[C]=CO1(25859)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC1C=[C]C([O])O1(25652)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[O]C1[CH]COC=C=C1(25753)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(94.1421,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 91.3 to 94.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC(=O)C=C=CO(25860)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC(O)C=C=C=O(25861)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC([O])=C[C]=C[O](25862)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=C([O])C[C]=C[O](25863)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=CC([O])[C]C=C[O](25864)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=CC([O])[CH]C=[C][O](25865)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C]C([O])C[C]=C[O](25866)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]CC([O])[C]=C=C[O](25867)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C=C([O])[CH]C=C[O](25868)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]C([O])[CH]C=C[O](25869)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C[CH]C([O])[C]=C=C[O](25870)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C=C([O])C=[C]C[O](25871)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=CC([O])C[C]=C[O](25872)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C]C([O])C=[C]C[O](25873)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=CC([O])[CH]C=C[O](25874)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]CC([O])[CH][C]=C=O(25875)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=CC([O])C=[C]C[O](25876)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C[CH]C([O])[CH][C]=C=O(25877)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([O])CC=C[O](22586)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.214e+09,'s^-1'), n=0.749, Ea=(200.242,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for hex_1_ene_5_yne
Exact match found for rate rule [hex_1_ene_5_yne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[O]C=C1[CH]C([O])[CH]C1(25878)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction52',
    reactants = ['O(T)(63)', '[CH2]C=CC=C=C[O](23648)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['O(T)(63)', '[CH]=C=CC([O])C=C(20057)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C1C(C=O)=CC1[O](25338)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_2H_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['H(8)', 'C=CC([O])C#CC=O(25879)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(39.9391,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 34.7 to 39.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=CC([O])[C]=CC=O(25880)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC([O])[CH]C=C=O(25881)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[CH2]C=C([O])C=CC=O(25882)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=[C]C([O])C=CC=O(25883)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]=CC([O])C=CC=O(25884)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['[O]C1[CH]CC(C=O)=C1(25670)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_pri_2H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC1[CH][C]=COO1(25602)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(122.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 118.5 to 122.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction64',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC(=O)C=CC=O(25885)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['C=CC([O])C=C=C[O](22584)'],
    products = ['C=CC1C=C(C=O)O1(25505)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]=O(373)', '[C]=CC([O])C=C(13789)'],
    products = ['C=CC([O])C=C=C[O](22584)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4885',
    isomers = [
        'C=CC([O])C=C=C[O](22584)',
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
    label = '4885',
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


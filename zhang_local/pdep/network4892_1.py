species(
    label = '[CH]=C=COO[CH]C=C(22591)',
    structure = SMILES('[CH]=C=COO[CH]C=C'),
    E0 = (442.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.48243,'amu*angstrom^2'), symmetry=1, barrier=(34.084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47957,'amu*angstrom^2'), symmetry=1, barrier=(34.0182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48138,'amu*angstrom^2'), symmetry=1, barrier=(34.0598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48032,'amu*angstrom^2'), symmetry=1, barrier=(34.0355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390171,0.0677683,-4.27675e-05,-4.12737e-10,6.91188e-12,53317.2,31.4896], Tmin=(100,'K'), Tmax=(986.066,'K')), NASAPolynomial(coeffs=[18.6527,0.0189194,-6.84395e-06,1.25159e-09,-8.97039e-14,48488.8,-62.5808], Tmin=(986.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=C=CJ)"""),
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
    label = '[CH]=[C]C1OOC1C=C(25579)',
    structure = SMILES('[CH]=[C]C1OOC1C=C'),
    E0 = (609.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930553,0.0514606,6.45993e-06,-5.91138e-08,3.12739e-11,73452.4,28.9834], Tmin=(100,'K'), Tmax=(900.682,'K')), NASAPolynomial(coeffs=[19.5229,0.0123784,-8.77507e-07,-7.51977e-11,7.00417e-15,68339.3,-68.5614], Tmin=(900.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1[CH]OOC=C=C1(25390)',
    structure = SMILES('[CH2]C1[CH]OOC=C=C1'),
    E0 = (540.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0693,0.0539737,-1.81841e-05,-1.33946e-08,8.42998e-12,65121.5,22.4124], Tmin=(100,'K'), Tmax=(1036.18,'K')), NASAPolynomial(coeffs=[13.5674,0.0273255,-1.08742e-05,2.01868e-09,-1.42251e-13,61371.9,-43.9231], Tmin=(1036.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCsJOOC) + radical(Isobutyl)"""),
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
    label = '[CH]=C=COOC=C=C(25580)',
    structure = SMILES('[CH]=C=COOC=C=C'),
    E0 = (498.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.32045,'amu*angstrom^2'), symmetry=1, barrier=(30.3598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3291,'amu*angstrom^2'), symmetry=1, barrier=(30.5586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31957,'amu*angstrom^2'), symmetry=1, barrier=(30.3396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15923,0.0776059,-8.20802e-05,4.16563e-08,-7.73838e-12,60133.6,30.1592], Tmin=(100,'K'), Tmax=(1015.77,'K')), NASAPolynomial(coeffs=[18.315,0.0159408,-5.53549e-06,9.4636e-10,-6.37915e-14,55938.1,-60.2124], Tmin=(1015.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=COOC[C]=C(25581)',
    structure = SMILES('[CH]=C=COOC[C]=C'),
    E0 = (562.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,350,500,795,815,540,610,2055,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(2.34529,'amu*angstrom^2'), symmetry=1, barrier=(53.9229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865093,'amu*angstrom^2'), symmetry=1, barrier=(19.8902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864759,'amu*angstrom^2'), symmetry=1, barrier=(19.8825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864484,'amu*angstrom^2'), symmetry=1, barrier=(19.8762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533072,0.0787694,-9.07028e-05,5.69015e-08,-1.44743e-11,67798,30.863], Tmin=(100,'K'), Tmax=(952.052,'K')), NASAPolynomial(coeffs=[12.2882,0.0293823,-1.28936e-05,2.41771e-09,-1.67736e-13,65559.6,-25.2712], Tmin=(952.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COOCC=[CH](25582)',
    structure = SMILES('[CH]=C=COOCC=[CH]'),
    E0 = (571.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,2750,2850,1437.5,1250,1305,750,350,540,610,2055,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.39539,'amu*angstrom^2'), symmetry=1, barrier=(32.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00333,'amu*angstrom^2'), symmetry=1, barrier=(23.0685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.996639,'amu*angstrom^2'), symmetry=1, barrier=(22.9147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.997899,'amu*angstrom^2'), symmetry=1, barrier=(22.9437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428511,0.0785116,-8.5721e-05,4.94318e-08,-1.14265e-11,68917,31.1628], Tmin=(100,'K'), Tmax=(1049.33,'K')), NASAPolynomial(coeffs=[14.1082,0.0263651,-1.11784e-05,2.07304e-09,-1.43403e-13,66046.1,-35.4927], Tmin=(1049.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C#COO[CH]C=C(25583)',
    structure = SMILES('[CH2]C#COO[CH]C=C'),
    E0 = (492.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,500,795,815,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2100,2250,500,550,458.205,458.293],'cm^-1')),
        HinderedRotor(inertia=(0.238498,'amu*angstrom^2'), symmetry=1, barrier=(35.5426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238475,'amu*angstrom^2'), symmetry=1, barrier=(35.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238489,'amu*angstrom^2'), symmetry=1, barrier=(35.5435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238495,'amu*angstrom^2'), symmetry=1, barrier=(35.5429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238421,'amu*angstrom^2'), symmetry=1, barrier=(35.5427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478223,0.0668925,-4.92488e-05,1.52233e-08,-1.20898e-12,59408.7,31.7872], Tmin=(100,'K'), Tmax=(1219.31,'K')), NASAPolynomial(coeffs=[17.3967,0.0233272,-1.0339e-05,1.97789e-09,-1.39414e-13,54395.6,-56.828], Tmin=(1219.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=CCJO) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]OOCC=C(25584)',
    structure = SMILES('[CH]=C=[C]OOCC=C'),
    E0 = (564.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,350,500,795,815,540,610,2055,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.464039,'amu*angstrom^2'), symmetry=1, barrier=(10.6692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.72144,'amu*angstrom^2'), symmetry=1, barrier=(85.5632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.464167,'amu*angstrom^2'), symmetry=1, barrier=(10.6721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.72773,'amu*angstrom^2'), symmetry=1, barrier=(85.7078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747063,0.0762452,-9.38267e-05,6.77166e-08,-2.03139e-11,68016.9,32.7514], Tmin=(100,'K'), Tmax=(805.12,'K')), NASAPolynomial(coeffs=[9.19497,0.0342726,-1.56253e-05,2.96061e-09,-2.05582e-13,66656.7,-6.17336], Tmin=(805.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C][CH]OOC=C=C(25585)',
    structure = SMILES('C=[C][CH]OOC=C=C'),
    E0 = (525.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32694,'amu*angstrom^2'), symmetry=1, barrier=(30.5091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32715,'amu*angstrom^2'), symmetry=1, barrier=(30.5139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32675,'amu*angstrom^2'), symmetry=1, barrier=(30.5045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32676,'amu*angstrom^2'), symmetry=1, barrier=(30.5049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169608,0.0731988,-6.48951e-05,2.81764e-08,-4.79616e-12,63350.5,32.143], Tmin=(100,'K'), Tmax=(1423.32,'K')), NASAPolynomial(coeffs=[19.4506,0.0190131,-7.79045e-06,1.42935e-09,-9.81769e-14,57861.9,-67.6826], Tmin=(1423.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]OOC=C=C(25586)',
    structure = SMILES('[CH]=C[CH]OOC=C=C'),
    E0 = (534.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.34309,'amu*angstrom^2'), symmetry=1, barrier=(30.8802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34348,'amu*angstrom^2'), symmetry=1, barrier=(30.8893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3414,'amu*angstrom^2'), symmetry=1, barrier=(30.8415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34127,'amu*angstrom^2'), symmetry=1, barrier=(30.8384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280291,0.0705153,-5.20075e-05,1.13761e-08,1.8063e-12,64460.1,31.6637], Tmin=(100,'K'), Tmax=(1054.15,'K')), NASAPolynomial(coeffs=[18.8746,0.0198469,-8.20868e-06,1.57409e-09,-1.13556e-13,59434.8,-64.265], Tmin=(1054.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=CCJO)"""),
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
    label = 'C=C[CH]O[O](6572)',
    structure = SMILES('C=C[CH]O[O]'),
    E0 = (193.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0263592,'amu*angstrom^2'), symmetry=1, barrier=(16.4069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264447,'amu*angstrom^2'), symmetry=1, barrier=(37.5506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44754,0.0263949,2.69871e-06,-2.29929e-08,1.07611e-11,23370.7,18.6111], Tmin=(100,'K'), Tmax=(985.371,'K')), NASAPolynomial(coeffs=[10.196,0.0131636,-4.89974e-06,9.15827e-10,-6.6401e-14,20959,-23.1458], Tmin=(985.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COO[CH][C]=C(25587)',
    structure = SMILES('[CH]=C=COO[CH][C]=C'),
    E0 = (679.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.48278,'amu*angstrom^2'), symmetry=1, barrier=(34.092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48344,'amu*angstrom^2'), symmetry=1, barrier=(34.1072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48312,'amu*angstrom^2'), symmetry=1, barrier=(34.0997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48326,'amu*angstrom^2'), symmetry=1, barrier=(34.1031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.186816,0.0742981,-7.39003e-05,3.64441e-08,-6.99214e-12,81928,32.6977], Tmin=(100,'K'), Tmax=(1278.54,'K')), NASAPolynomial(coeffs=[18.6835,0.01643,-6.00863e-06,1.04353e-09,-7.00684e-14,77198.2,-61.083], Tmin=(1278.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=COO[CH]C=[CH](25588)',
    structure = SMILES('[CH]=C=COO[CH]C=[CH]'),
    E0 = (689.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,3025,407.5,1350,352.5,540,610,2055,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.48663,'amu*angstrom^2'), symmetry=1, barrier=(34.1806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48763,'amu*angstrom^2'), symmetry=1, barrier=(34.2035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48619,'amu*angstrom^2'), symmetry=1, barrier=(34.1705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48547,'amu*angstrom^2'), symmetry=1, barrier=(34.154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311843,0.0713463,-5.96248e-05,1.72377e-08,8.77955e-13,83037.1,32.1738], Tmin=(100,'K'), Tmax=(992.527,'K')), NASAPolynomial(coeffs=[18.9141,0.0160504,-5.78901e-06,1.04787e-09,-7.44248e-14,78375.4,-62.3135], Tmin=(992.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=C=CJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C=[C]OO[CH]C=C(25589)',
    structure = SMILES('[CH]=C=[C]OO[CH]C=C'),
    E0 = (681.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,540,610,2055,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,352.099],'cm^-1')),
        HinderedRotor(inertia=(0.289801,'amu*angstrom^2'), symmetry=1, barrier=(25.9305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.89983,'amu*angstrom^2'), symmetry=1, barrier=(79.2125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292468,'amu*angstrom^2'), symmetry=1, barrier=(25.7754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901459,'amu*angstrom^2'), symmetry=1, barrier=(79.288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.710142,0.0681908,-6.47412e-05,3.15346e-08,-6.13601e-12,82133.3,33.4716], Tmin=(100,'K'), Tmax=(1238.41,'K')), NASAPolynomial(coeffs=[14.8614,0.0224829,-9.37844e-06,1.73146e-09,-1.19585e-13,78628.3,-37.826], Tmin=(1238.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=CCJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COOC1[CH]C1(25590)',
    structure = SMILES('[CH]=C=COOC1[CH]C1'),
    E0 = (541.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491374,0.0658129,-4.01184e-05,-2.12646e-09,7.4098e-12,65308.5,31.1845], Tmin=(100,'K'), Tmax=(980.604,'K')), NASAPolynomial(coeffs=[18.1388,0.018892,-6.68591e-06,1.21064e-09,-8.64435e-14,60642.4,-59.753], Tmin=(980.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(C=C=CJ) + radical(CCJCOOH)"""),
)

species(
    label = 'C=C[CH]OOC1[C]=C1(25591)',
    structure = SMILES('C=C[CH]OOC1[C]=C1'),
    E0 = (585.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.566785,0.0632134,-3.48615e-05,-2.10146e-09,5.23778e-12,70586.2,30.4893], Tmin=(100,'K'), Tmax=(1086.09,'K')), NASAPolynomial(coeffs=[17.7775,0.0225718,-1.01435e-05,2.00744e-09,-1.46272e-13,65506.3,-60.1394], Tmin=(1086.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CCJO) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = '[CH]C1=COOC1C=C(25592)',
    structure = SMILES('[CH]C1=COOC1C=C'),
    E0 = (395.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927274,0.050409,1.38833e-05,-5.58136e-08,2.50184e-11,47658.6,28.6864], Tmin=(100,'K'), Tmax=(987.751,'K')), NASAPolynomial(coeffs=[17.1702,0.0253193,-9.80396e-06,1.8768e-09,-1.38162e-13,42465,-59.5233], Tmin=(987.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C1=CC[CH][CH]OOC=1(25593)',
    structure = SMILES('C1=CC[CH][CH]OOC=1'),
    E0 = (508.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80419,0.0304462,5.06067e-05,-8.29643e-08,3.20905e-11,61261.5,23.9613], Tmin=(100,'K'), Tmax=(1000.54,'K')), NASAPolynomial(coeffs=[13.4392,0.0278692,-1.14014e-05,2.24278e-09,-1.66546e-13,56733.9,-43.1678], Tmin=(1000.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclooctane) + radical(CCJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = 'C=C=COOC=C=C(25594)',
    structure = SMILES('C=C=COOC=C=C'),
    E0 = (344.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0998726,0.0770565,-7.52051e-05,3.63685e-08,-6.87533e-12,41557.8,29.7516], Tmin=(100,'K'), Tmax=(1291.09,'K')), NASAPolynomial(coeffs=[18.9466,0.0186654,-7.36446e-06,1.33764e-09,-9.20114e-14,36691.4,-65.9874], Tmin=(1291.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C1=COC=1(22275)',
    structure = SMILES('C1=COC=1'),
    E0 = (388.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17882,0.00883829,2.87655e-05,-4.70958e-08,1.95464e-11,46730.5,10.2865], Tmin=(100,'K'), Tmax=(941.322,'K')), NASAPolynomial(coeffs=[9.87744,0.00380291,-5.45265e-07,1.04211e-10,-1.15627e-14,44431.4,-27.139], Tmin=(941.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH]=C(C=O)O[CH]C=C(22593)',
    structure = SMILES('[CH]=C(C=O)O[CH]C=C'),
    E0 = (160.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,373.034,373.043,373.044],'cm^-1')),
        HinderedRotor(inertia=(0.206408,'amu*angstrom^2'), symmetry=1, barrier=(20.3848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206417,'amu*angstrom^2'), symmetry=1, barrier=(20.3849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206432,'amu*angstrom^2'), symmetry=1, barrier=(20.3847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20643,'amu*angstrom^2'), symmetry=1, barrier=(20.3847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3744.15,'J/mol'), sigma=(6.10363,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.83 K, Pc=37.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611721,0.0731945,-6.82978e-05,3.16834e-08,-5.90605e-12,19466.9,27.1732], Tmin=(100,'K'), Tmax=(1275.67,'K')), NASAPolynomial(coeffs=[15.7913,0.0255974,-1.23307e-05,2.43494e-09,-1.74084e-13,15594.1,-49.7552], Tmin=(1275.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=CO[O](20803)',
    structure = SMILES('[CH]=C=CO[O]'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,180,180,180,2867.84],'cm^-1')),
        HinderedRotor(inertia=(0.0105329,'amu*angstrom^2'), symmetry=1, barrier=(9.09612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13751,0.0383774,-4.88548e-05,3.09451e-08,-7.37348e-12,48769,18.1869], Tmin=(100,'K'), Tmax=(1167.31,'K')), NASAPolynomial(coeffs=[10.2545,0.00577971,-8.1986e-07,1.20449e-12,5.54427e-15,47199.8,-20.8326], Tmin=(1167.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]C=C(18735)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,192.655,193.544,193.915],'cm^-1')),
        HinderedRotor(inertia=(1.88068,'amu*angstrom^2'), symmetry=1, barrier=(50.3487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32096,0.00806329,3.46645e-05,-4.52343e-08,1.64854e-11,45350.1,10.7121], Tmin=(100,'K'), Tmax=(975.253,'K')), NASAPolynomial(coeffs=[5.21066,0.0176207,-6.65616e-06,1.20944e-09,-8.49962e-14,44158.4,-2.57721], Tmin=(975.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]OOC=C=[CH](23157)',
    structure = SMILES('[CH]OOC=C=[CH]'),
    E0 = (713.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,3010,987.5,1337.5,450,1655,311.097,851.355],'cm^-1')),
        HinderedRotor(inertia=(0.594343,'amu*angstrom^2'), symmetry=1, barrier=(13.6651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59463,'amu*angstrom^2'), symmetry=1, barrier=(13.6717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.84824,'amu*angstrom^2'), symmetry=1, barrier=(88.4785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31984,0.0585742,-7.51673e-05,4.76862e-08,-1.1737e-11,85949.5,22.4565], Tmin=(100,'K'), Tmax=(1002.18,'K')), NASAPolynomial(coeffs=[13.0987,0.0115612,-4.80092e-06,8.7714e-10,-6.01566e-14,83588.6,-34.3952], Tmin=(1002.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CH2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COO[C]=C[CH2](25595)',
    structure = SMILES('[CH]=C=COO[C]=C[CH2]'),
    E0 = (713.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.939183,'amu*angstrom^2'), symmetry=1, barrier=(21.5937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939321,'amu*angstrom^2'), symmetry=1, barrier=(21.5968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43883,'amu*angstrom^2'), symmetry=1, barrier=(56.0736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938624,'amu*angstrom^2'), symmetry=1, barrier=(21.5808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510395,0.0750457,-8.35094e-05,4.81156e-08,-1.09678e-11,85931.6,33.2578], Tmin=(100,'K'), Tmax=(1071.25,'K')), NASAPolynomial(coeffs=[14.8126,0.0216423,-8.73296e-06,1.5807e-09,-1.07992e-13,82867.3,-36.7267], Tmin=(1071.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[C]#C[CH]OO[CH]C=C(25596)',
    structure = SMILES('[C]#C[CH]OO[CH]C=C'),
    E0 = (803.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,350,500,795,815,2175,525,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742431,0.0722311,-7.6891e-05,4.44442e-08,-1.04341e-11,96726.5,29.8198], Tmin=(100,'K'), Tmax=(1027.83,'K')), NASAPolynomial(coeffs=[12.2222,0.0275556,-1.16925e-05,2.15564e-09,-1.4829e-13,94366.6,-25.8785], Tmin=(1027.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl) + radical(C=CCJO)"""),
)

species(
    label = 'C#CC1OO[CH]C1[CH2](25477)',
    structure = SMILES('C#CC1OO[CH]C1[CH2]'),
    E0 = (467.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04411,0.0491009,1.79811e-05,-7.80688e-08,4.16265e-11,56345.9,25.3711], Tmin=(100,'K'), Tmax=(856.128,'K')), NASAPolynomial(coeffs=[19.6393,0.010011,2.73672e-06,-9.95327e-10,8.00809e-14,51410.5,-71.6803], Tmin=(856.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxolane) + radical(CCsJOOC) + radical(Isobutyl)"""),
)

species(
    label = '[C]#CCOO[CH]C=C(25597)',
    structure = SMILES('[C]#CCOO[CH]C=C'),
    E0 = (616.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2175,525,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.610027,0.0719163,-6.83394e-05,3.47061e-08,-7.14367e-12,74310.3,29.9952], Tmin=(100,'K'), Tmax=(1166.34,'K')), NASAPolynomial(coeffs=[13.5742,0.027455,-1.11587e-05,2.02207e-09,-1.37965e-13,71286.1,-34.5443], Tmin=(1166.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(616.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CCJO) + radical(Acetyl)"""),
)

species(
    label = 'C#CCOO[CH][C]=C(25598)',
    structure = SMILES('C#CCOO[CH][C]=C'),
    E0 = (517.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,1685,370,2175,525,180,435.869],'cm^-1')),
        HinderedRotor(inertia=(0.551768,'amu*angstrom^2'), symmetry=1, barrier=(74.3869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213785,'amu*angstrom^2'), symmetry=1, barrier=(28.8205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551777,'amu*angstrom^2'), symmetry=1, barrier=(74.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551778,'amu*angstrom^2'), symmetry=1, barrier=(74.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213777,'amu*angstrom^2'), symmetry=1, barrier=(28.8204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619078,0.0721989,-6.77846e-05,3.33717e-08,-6.65847e-12,62366.2,29.6947], Tmin=(100,'K'), Tmax=(1198.1,'K')), NASAPolynomial(coeffs=[14.0255,0.0274395,-1.17461e-05,2.18959e-09,-1.51834e-13,59153.8,-37.4067], Tmin=(1198.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C[CH]OOCC#C(25599)',
    structure = SMILES('[CH]=C[CH]OOCC#C'),
    E0 = (526.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,3025,407.5,1350,352.5,350,500,795,815,2175,525,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,351.45],'cm^-1')),
        HinderedRotor(inertia=(0.375804,'amu*angstrom^2'), symmetry=1, barrier=(32.9371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714654,'amu*angstrom^2'), symmetry=1, barrier=(62.6374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375782,'amu*angstrom^2'), symmetry=1, barrier=(32.9373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375812,'amu*angstrom^2'), symmetry=1, barrier=(32.937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714606,'amu*angstrom^2'), symmetry=1, barrier=(62.6371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367593,0.0735778,-6.80875e-05,3.21711e-08,-6.05279e-12,63491.8,30.528], Tmin=(100,'K'), Tmax=(1282.16,'K')), NASAPolynomial(coeffs=[16.4044,0.0235477,-9.55793e-06,1.73861e-09,-1.19023e-13,59379.4,-50.8264], Tmin=(1282.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=CCJO)"""),
)

species(
    label = '[C]#C[CH]OOCC=C(25600)',
    structure = SMILES('[C]#C[CH]OOCC=C'),
    E0 = (685.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2175,525,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474539,0.0840738,-0.000120079,1.0003e-07,-3.33526e-11,82622.9,30.1794], Tmin=(100,'K'), Tmax=(849.342,'K')), NASAPolynomial(coeffs=[8.16868,0.036586,-1.63401e-05,3.00546e-09,-2.02722e-13,81721.8,-3.29483], Tmin=(849.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl)"""),
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
    label = 'C#CCOOC=C=C(25603)',
    structure = SMILES('C#CCOOC=C=C'),
    E0 = (336.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.422177,0.0774002,-8.21202e-05,4.59353e-08,-1.02824e-11,40579.4,27.7706], Tmin=(100,'K'), Tmax=(1083.75,'K')), NASAPolynomial(coeffs=[14.4221,0.0257274,-1.05995e-05,1.93879e-09,-1.33152e-13,37545,-40.8966], Tmin=(1083.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1OOC1C=C(22598)',
    structure = SMILES('C#CC1OOC1C=C'),
    E0 = (290.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0167,0.0461175,2.95624e-05,-9.02524e-08,4.48873e-11,35087.3,25.417], Tmin=(100,'K'), Tmax=(884.326,'K')), NASAPolynomial(coeffs=[21.4339,0.00785252,2.72579e-06,-8.59502e-10,6.37691e-14,29361.4,-82.5318], Tmin=(884.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane)"""),
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
    label = '[CH]OO[CH]C=C(13499)',
    structure = SMILES('[CH]OO[CH]C=C'),
    E0 = (502.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,500,795,815,3010,987.5,1337.5,450,1655,180,180,1029.59,1029.98],'cm^-1')),
        HinderedRotor(inertia=(0.00288069,'amu*angstrom^2'), symmetry=1, barrier=(2.16641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2133,'amu*angstrom^2'), symmetry=1, barrier=(27.8963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370475,'amu*angstrom^2'), symmetry=1, barrier=(27.8928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21299,'amu*angstrom^2'), symmetry=1, barrier=(27.8891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35175,0.049849,-3.4844e-05,7.97716e-09,5.10098e-13,60563.3,23.8792], Tmin=(100,'K'), Tmax=(1149.79,'K')), NASAPolynomial(coeffs=[14.4637,0.0165748,-7.5343e-06,1.47719e-09,-1.06227e-13,56732.4,-44.7556], Tmin=(1149.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C=COOC=[C]C(25604)',
    structure = SMILES('[CH]=C=COOC=[C]C'),
    E0 = (560.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261313,0.0804718,-9.03753e-05,5.26808e-08,-1.21374e-11,67490.7,31.2159], Tmin=(100,'K'), Tmax=(1061.02,'K')), NASAPolynomial(coeffs=[15.5017,0.0230161,-9.14803e-06,1.64339e-09,-1.11775e-13,64256.6,-43.2126], Tmin=(1061.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=COO[C]=CC(25605)',
    structure = SMILES('[CH]=C=COO[C]=CC'),
    E0 = (561.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635942,0.0759852,-8.63114e-05,5.36672e-08,-1.35309e-11,67702.8,32.5327], Tmin=(100,'K'), Tmax=(960.695,'K')), NASAPolynomial(coeffs=[12.0278,0.0285532,-1.22521e-05,2.27407e-09,-1.56889e-13,65514,-21.9696], Tmin=(960.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=[C]OOC=CC(25606)',
    structure = SMILES('[CH]=C=[C]OOC=CC'),
    E0 = (561.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,500,795,815,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.793457,'amu*angstrom^2'), symmetry=1, barrier=(18.2431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793293,'amu*angstrom^2'), symmetry=1, barrier=(18.2394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793481,'amu*angstrom^2'), symmetry=1, barrier=(18.2437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793213,'amu*angstrom^2'), symmetry=1, barrier=(18.2375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635942,0.0759852,-8.63114e-05,5.36672e-08,-1.35309e-11,67702.8,32.5327], Tmin=(100,'K'), Tmax=(960.695,'K')), NASAPolynomial(coeffs=[12.0278,0.0285532,-1.22521e-05,2.27407e-09,-1.56889e-13,65514,-21.9696], Tmin=(960.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]OOC=C=C(25607)',
    structure = SMILES('[CH2]C=[C]OOC=C=C'),
    E0 = (558.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,350,500,795,815,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.858718,'amu*angstrom^2'), symmetry=1, barrier=(19.7436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861388,'amu*angstrom^2'), symmetry=1, barrier=(19.805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41859,'amu*angstrom^2'), symmetry=1, barrier=(32.6162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85891,'amu*angstrom^2'), symmetry=1, barrier=(19.748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625523,0.0725315,-7.01985e-05,3.50986e-08,-7.06926e-12,67348.1,32.2181], Tmin=(100,'K'), Tmax=(1190.31,'K')), NASAPolynomial(coeffs=[14.5122,0.0258657,-1.13914e-05,2.16196e-09,-1.51601e-13,64042.2,-37.1964], Tmin=(1190.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1[CH]OOC=CC1(25608)',
    structure = SMILES('[CH]=C1[CH]OOC=CC1'),
    E0 = (426.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45383,-0.00773859,0.000257,-3.77923e-07,1.59194e-10,51457.5,27.5308], Tmin=(100,'K'), Tmax=(906.019,'K')), NASAPolynomial(coeffs=[45.3542,-0.034955,2.62368e-05,-5.16799e-09,3.36641e-13,36664.8,-217.666], Tmin=(906.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]COO[C]=C[CH2](25609)',
    structure = SMILES('[CH]=[C]COO[C]=C[CH2]'),
    E0 = (869.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,289.139],'cm^-1')),
        HinderedRotor(inertia=(0.105202,'amu*angstrom^2'), symmetry=1, barrier=(6.24107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105202,'amu*angstrom^2'), symmetry=1, barrier=(6.24106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105202,'amu*angstrom^2'), symmetry=1, barrier=(6.24107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72314,'amu*angstrom^2'), symmetry=1, barrier=(42.9008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.723166,'amu*angstrom^2'), symmetry=1, barrier=(42.9007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576646,0.0814637,-0.000111153,8.97597e-08,-3.0148e-11,104747,34.8347], Tmin=(100,'K'), Tmax=(736.852,'K')), NASAPolynomial(coeffs=[8.65712,0.0366006,-1.7794e-05,3.45432e-09,-2.42404e-13,103584,-1.49753], Tmin=(736.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(869.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]OOC[CH][CH2](25610)',
    structure = SMILES('[CH]=C=[C]OOC[CH][CH2]'),
    E0 = (836.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.274944,'amu*angstrom^2'), symmetry=1, barrier=(6.32151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27496,'amu*angstrom^2'), symmetry=1, barrier=(6.32188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274963,'amu*angstrom^2'), symmetry=1, barrier=(6.32193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0418649,'amu*angstrom^2'), symmetry=1, barrier=(57.6907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.50917,'amu*angstrom^2'), symmetry=1, barrier=(57.6907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25176,0.0900607,-0.000137218,1.16963e-07,-3.94549e-11,100785,37.7522], Tmin=(100,'K'), Tmax=(840.419,'K')), NASAPolynomial(coeffs=[9.12194,0.0351013,-1.63838e-05,3.07066e-09,-2.0906e-13,99743.5,-0.822129], Tmin=(840.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO) + radical(RCCJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=C[CH]OO[C]=C[CH2](25611)',
    structure = SMILES('[CH]=C[CH]OO[C]=C[CH2]'),
    E0 = (749.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,500,795,815,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,450.959],'cm^-1')),
        HinderedRotor(inertia=(0.221041,'amu*angstrom^2'), symmetry=1, barrier=(31.9356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0163258,'amu*angstrom^2'), symmetry=1, barrier=(2.35651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22111,'amu*angstrom^2'), symmetry=1, barrier=(31.9363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221147,'amu*angstrom^2'), symmetry=1, barrier=(31.9358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464061,0.0698572,-5.96754e-05,2.52146e-08,-4.23109e-12,90265.3,35.3665], Tmin=(100,'K'), Tmax=(1425.9,'K')), NASAPolynomial(coeffs=[17.4437,0.0222251,-9.56788e-06,1.78726e-09,-1.23617e-13,85423,-52.5749], Tmin=(1425.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(749.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=[C]OO[CH]C[CH2](25612)',
    structure = SMILES('[CH]=C=[C]OO[CH]C[CH2]'),
    E0 = (822.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,196.305],'cm^-1')),
        HinderedRotor(inertia=(0.254016,'amu*angstrom^2'), symmetry=1, barrier=(6.95066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00591,'amu*angstrom^2'), symmetry=1, barrier=(54.8696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253967,'amu*angstrom^2'), symmetry=1, barrier=(6.95048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00842,'amu*angstrom^2'), symmetry=1, barrier=(54.8727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133218,'amu*angstrom^2'), symmetry=1, barrier=(54.8712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258045,0.089592,-0.000134406,1.13616e-07,-3.81777e-11,99103.9,35.7201], Tmin=(100,'K'), Tmax=(837.027,'K')), NASAPolynomial(coeffs=[9.16908,0.0353883,-1.64468e-05,3.08024e-09,-2.09825e-13,98019.2,-3.25371], Tmin=(837.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(822.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=CJO) + radical(C=C=CJ) + radical(CCsJOOC)"""),
)

species(
    label = '[CH]=[C]COO[CH][C]=C(25613)',
    structure = SMILES('[CH]=[C]COO[CH][C]=C'),
    E0 = (836.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440,180,609.49],'cm^-1')),
        HinderedRotor(inertia=(0.115444,'amu*angstrom^2'), symmetry=1, barrier=(2.65429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34407,'amu*angstrom^2'), symmetry=1, barrier=(30.9029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163659,'amu*angstrom^2'), symmetry=1, barrier=(43.3899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0536724,'amu*angstrom^2'), symmetry=1, barrier=(14.0221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88758,'amu*angstrom^2'), symmetry=1, barrier=(43.3993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904979,0.0729217,-7.36685e-05,4.06959e-08,-9.45029e-12,100716,31.9425], Tmin=(100,'K'), Tmax=(1014.23,'K')), NASAPolynomial(coeffs=[10.6534,0.034475,-1.68075e-05,3.32022e-09,-2.37445e-13,98738.1,-15.2257], Tmin=(1014.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]OO[CH][C]=C(25614)',
    structure = SMILES('[CH]=C[CH]OO[CH][C]=C'),
    E0 = (715.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3000,3050,390,425,1340,1360,335,370,350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655,391.729,391.737],'cm^-1')),
        HinderedRotor(inertia=(0.361952,'amu*angstrom^2'), symmetry=1, barrier=(39.4146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361951,'amu*angstrom^2'), symmetry=1, barrier=(39.4145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361962,'amu*angstrom^2'), symmetry=1, barrier=(39.4145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00109857,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361955,'amu*angstrom^2'), symmetry=1, barrier=(39.4144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459241,0.0654815,-3.80902e-05,-8.09865e-10,5.31422e-12,86247.8,33.6555], Tmin=(100,'K'), Tmax=(1065.05,'K')), NASAPolynomial(coeffs=[18.2524,0.0219478,-9.58201e-06,1.87893e-09,-1.36752e-13,81136.6,-59.5092], Tmin=(1065.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(715.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CCJO) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C1=COC1(5259)',
    structure = SMILES('C1=COC1'),
    E0 = (4.81952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10338,-0.00295546,9.94696e-05,-1.38902e-07,5.64804e-11,633.036,9.26427], Tmin=(100,'K'), Tmax=(918.619,'K')), NASAPolynomial(coeffs=[16.7387,-0.00374949,5.11332e-06,-1.00759e-09,6.07118e-14,-4343.72,-68.8138], Tmin=(918.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.81952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C=COC([CH2])C=O(22595)',
    structure = SMILES('[CH]=C=COC([CH2])C=O'),
    E0 = (203.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06237,'amu*angstrom^2'), symmetry=1, barrier=(24.4259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06109,'amu*angstrom^2'), symmetry=1, barrier=(24.3965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06115,'amu*angstrom^2'), symmetry=1, barrier=(24.3979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06101,'amu*angstrom^2'), symmetry=1, barrier=(24.3948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3838.73,'J/mol'), sigma=(6.28681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=599.60 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11754,0.0888871,-9.45669e-05,4.71822e-08,-8.70488e-12,24670.9,34.6888], Tmin=(100,'K'), Tmax=(1541.19,'K')), NASAPolynomial(coeffs=[25.0164,0.00605141,6.61496e-07,-3.28778e-10,2.68985e-14,18397.8,-96.9145], Tmin=(1541.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C1[CH]OO[CH]C1[CH2](25365)',
    structure = SMILES('[CH]=C1[CH]OO[CH]C1[CH2]'),
    E0 = (651.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09545,0.0411671,4.39336e-05,-8.95045e-08,3.69232e-11,78526.8,25.7999], Tmin=(100,'K'), Tmax=(993.662,'K')), NASAPolynomial(coeffs=[19.8837,0.0208102,-8.77846e-06,1.84404e-09,-1.44666e-13,72064.1,-78.4545], Tmin=(993.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(Cds_P) + radical(C=CCJO) + radical(CCsJOOC) + radical(Isobutyl)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=COOC=[CH](23820)',
    structure = SMILES('[CH]=C=COOC=[CH]'),
    E0 = (605.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,540,610,2055,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.27599,'amu*angstrom^2'), symmetry=1, barrier=(29.3375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2661,'amu*angstrom^2'), symmetry=1, barrier=(29.1101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506102,0.0681032,-7.67651e-05,4.20328e-08,-8.7654e-12,72933.7,28.4577], Tmin=(100,'K'), Tmax=(1250.38,'K')), NASAPolynomial(coeffs=[17.9426,0.00916525,-2.27258e-06,2.95591e-10,-1.66214e-14,68820.1,-58.5723], Tmin=(1250.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CCOO[C]=C[CH2](25615)',
    structure = SMILES('C#CCOO[C]=C[CH2]'),
    E0 = (550.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,500,795,815,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2175,525,1104.7],'cm^-1')),
        HinderedRotor(inertia=(0.823848,'amu*angstrom^2'), symmetry=1, barrier=(18.9419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87323,'amu*angstrom^2'), symmetry=1, barrier=(66.0613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386384,'amu*angstrom^2'), symmetry=1, barrier=(8.88373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87627,'amu*angstrom^2'), symmetry=1, barrier=(66.1311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87362,'amu*angstrom^2'), symmetry=1, barrier=(66.0702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711555,0.0754951,-8.55356e-05,5.46356e-08,-1.43643e-11,66380.2,31.0961], Tmin=(100,'K'), Tmax=(916.809,'K')), NASAPolynomial(coeffs=[10.7304,0.0317822,-1.40146e-05,2.62706e-09,-1.81989e-13,64543.2,-16.3687], Tmin=(916.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[C]#C[CH]OOC=CC(25616)',
    structure = SMILES('[C]#C[CH]OOC=CC'),
    E0 = (683.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,350,500,795,815,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,273.374,273.424],'cm^-1')),
        HinderedRotor(inertia=(1.20521,'amu*angstrom^2'), symmetry=1, barrier=(63.9602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20503,'amu*angstrom^2'), symmetry=1, barrier=(63.9601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162302,'amu*angstrom^2'), symmetry=1, barrier=(8.61504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162375,'amu*angstrom^2'), symmetry=1, barrier=(8.61643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20521,'amu*angstrom^2'), symmetry=1, barrier=(63.9602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371727,0.0836464,-0.000111652,8.43564e-08,-2.56935e-11,82308.5,29.9353], Tmin=(100,'K'), Tmax=(842.568,'K')), NASAPolynomial(coeffs=[10.8023,0.0312339,-1.31913e-05,2.3742e-09,-1.5876e-13,80653.5,-17.99], Tmin=(842.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCsJOOC)"""),
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
    E0 = (442.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (609.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (567.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (723.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (442.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (715.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (703.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (658.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (608.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (770.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (725.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (442.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (700.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (891.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (901.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (893.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (668.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (592.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (598.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (508.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (467.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (525.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (755.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (786.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1037.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (925.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1015.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (550.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (769.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (641.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (571.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (776.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (505.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (517.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (481.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (450.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (1088.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (649.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (697.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (642.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (804.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (502.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (909.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (876.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (788.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (862.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (861.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (740.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (532.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (755.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (651.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (1020.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (442.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (595.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (811.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (449.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C=CC=O(5269)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=[C]C1OOC1C=C(25579)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(413769,'s^-1'), n=1.87624, Ea=(167.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 165.5 to 167.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH2]C1[CH]OOC=C=C1(25390)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra_cdsingleH] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=COOC=C=C(25580)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1092.27,'m^3/(mol*s)'), n=1.64867, Ea=(13.1815,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CC=O(5269)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(220.204,'m^3/(mol*s)'), n=1.265, Ea=(253.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_rad/OneDe]
Euclidian distance = 4.0
family: R_Addition_MultipleBond
Ea raised from 248.7 to 253.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=COOC[C]=C(25581)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.09427e+10,'s^-1'), n=1.04582, Ea=(152.506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=COOCC=[CH](25582)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.65416e+07,'s^-1'), n=1.654, Ea=(131.736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C#COO[CH]C=C(25583)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C=[C]OOCC=C(25584)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C][CH]OOC=C=C(25585)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R7HJ_1;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C[CH]OOC=C=C(25586)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R8Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=C[O](5266)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(81.9698,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 82.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[CH]O[O](6572)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.15767e+07,'m^3/(mol*s)'), n=0.0716491, Ea=(15.4197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_allenic;O_rad/NonDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', '[CH]=C=COO[CH][C]=C(25587)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH]=C=COO[CH]C=[CH](25588)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]=C=[C]OO[CH]C=C(25589)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C=COOC1[CH]C1(25590)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_csHO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C=C[CH]OOC1[C]=C1(25591)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(149.998,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]C1=COOC1C=C(25592)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C1=CC[CH][CH]OOC=1(25593)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(66.4311,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra_cdsingleH] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 61.4 to 66.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C=C=COOC=C=C(25594)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH2]C=C[O](5266)', 'C1=COC=1(22275)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO;Cd_pri_rad_in;OO_intra]
Euclidian distance = 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C(C=O)O[CH]C=C(22593)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C=CO[O](20803)', '[CH]C=C(18735)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(64)', '[CH]OOC=C=[CH](23157)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]=C=COO[C]=C[CH2](25595)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', '[C]#C[CH]OO[CH]C=C(25596)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C#CC1OO[CH]C1[CH2](25477)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;doublebond_intra_2H_pri;radadd_intra_csHDe] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[C]#CCOO[CH]C=C(25597)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.9263e+09,'s^-1'), n=1.08337, Ea=(153.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CCOO[CH][C]=C(25598)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.04875e+10,'s^-1'), n=0.94, Ea=(123.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C[CH]OOCC#C(25599)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]#C[CH]OOCC=C(25600)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(134600,'s^-1'), n=1.95079, Ea=(90.5935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C#CC1C[CH][CH]OO1(25601)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.0689e+08,'s^-1'), n=0.637531, Ea=(63.1312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHDe] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHCt]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C=CC1[CH][C]=COO1(25602)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(75.7068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCd] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C#CCOOC=C=C(25603)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C#CC1OOC1C=C(22598)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]#C(5143)', '[CH]OO[CH]C=C(13499)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C=COOC=[C]C(25604)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C=COO[C]=CC(25605)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C=[C]OOC=CC(25606)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(37753.8,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C=[C]OOC=C=C(25607)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C1[CH]OOC=CC1(25608)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(523137,'s^-1'), n=1.06759, Ea=(60.8507,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]COO[C]=C[CH2](25609)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C=[C]OOC[CH][CH2](25610)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C[CH]OO[C]=C[CH2](25611)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C=[C]OO[CH]C[CH2](25612)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]COO[CH][C]=C(25613)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C[CH]OO[CH][C]=C(25614)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C1=COC1(5259)', '[CH]=C=C[O](8556)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1','*|/',1.74), n=0, Ea=(89.9403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OOR] for rate rule [R3OO_SD;C_pri_rad_intra;OOR]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['[CH]=C=COC([CH2])C=O(22595)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C1[CH]OO[CH]C1[CH2](25365)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CH2(T)(28)', '[CH]=C=COOC=[CH](23820)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C=C[O](5266)', 'C#CC=O(21959)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(372.961,'m^3/(mol*s)'), n=1.215, Ea=(267.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Od_CO-CtH;YJ] for rate rule [Od_CO-CtH;O_rad/OneDe]
Euclidian distance = 4.0
family: R_Addition_MultipleBond
Ea raised from 262.7 to 267.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['C#CCOO[C]=C[CH2](25615)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/Ct]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[C]#C[CH]OOC=CC(25616)'],
    products = ['[CH]=C=COO[CH]C=C(22591)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.87683e+06,'s^-1'), n=1.52702, Ea=(127.963,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;Cs_H_out_2H] + [R8Hall;Y_rad_out;Cs_H_out] for rate rule [R8Hall;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=C=COO[CH]C=C(22591)'],
    products = ['C#CC1CC=COO1(25617)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""From training reaction 1 used for R6_SSSDS;C_rad_out_H/OneDe;Cpri_rad_out_2H
Exact match found for rate rule [R6_SSSDS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

network(
    label = '4892',
    isomers = [
        '[CH]=C=COO[CH]C=C(22591)',
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
    label = '4892',
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


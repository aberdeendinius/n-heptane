species(
    label = '[O]C=C=CC[C]=CO(25689)',
    structure = SMILES('[O]C=C=CC[C]=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CO(12571)',
    structure = SMILES('C=C=CO'),
    E0 = (-26.0646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.34368,'amu*angstrom^2'), symmetry=1, barrier=(30.8938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3437.21,'J/mol'), sigma=(5.57865,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.88 K, Pc=44.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31583,0.0236137,2.05754e-05,-5.73733e-08,2.79863e-11,-3061.58,12.125], Tmin=(100,'K'), Tmax=(901.949,'K')), NASAPolynomial(coeffs=[16.2977,-0.00239911,3.975e-06,-8.57293e-10,5.72973e-14,-7047.88,-62.0029], Tmin=(901.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.0646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[O]C=[C]C1CC1=CO(29139)',
    structure = SMILES('[O]C=[C]C1CC1=CO'),
    E0 = (220.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0215818,0.0685472,-2.21542e-05,-4.33127e-08,2.83969e-11,26664.8,25.6982], Tmin=(100,'K'), Tmax=(912.308,'K')), NASAPolynomial(coeffs=[26.5571,0.002819,2.69274e-06,-6.55026e-10,4.24508e-14,19716.7,-111.429], Tmin=(912.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = '[O]C=C=CC=C=CO(29140)',
    structure = SMILES('[O]C=C=CC=C=CO'),
    E0 = (98.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3671,'amu*angstrom^2'), symmetry=1, barrier=(31.4324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37376,'amu*angstrom^2'), symmetry=1, barrier=(31.5854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.21745,0.0721515,-2.80775e-05,-4.45686e-08,3.0974e-11,12072.9,25.9829], Tmin=(100,'K'), Tmax=(904.916,'K')), NASAPolynomial(coeffs=[29.885,-0.00475961,6.33532e-06,-1.35059e-09,9.03092e-14,4325.82,-128.939], Tmin=(904.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CCC#CO(29141)',
    structure = SMILES('[O]C=C=CCC#CO'),
    E0 = (189.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2100,2250,500,550,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13735,'amu*angstrom^2'), symmetry=1, barrier=(26.15,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13814,'amu*angstrom^2'), symmetry=1, barrier=(26.168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13816,'amu*angstrom^2'), symmetry=1, barrier=(26.1685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60722,0.0649476,-5.01777e-05,1.16249e-08,1.93342e-12,22971.9,27.6345], Tmin=(100,'K'), Tmax=(1008.66,'K')), NASAPolynomial(coeffs=[17.6274,0.0163932,-6.14014e-06,1.13627e-09,-8.13596e-14,18574.8,-59.4012], Tmin=(1008.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=COJ)"""),
)

species(
    label = 'O=C=C=CC[C]=CO(29142)',
    structure = SMILES('O=C=C=CC[C]=CO'),
    E0 = (235.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.124449,'amu*angstrom^2'), symmetry=1, barrier=(2.86133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122383,'amu*angstrom^2'), symmetry=1, barrier=(2.81383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.8502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20877,0.0489246,-1.80401e-05,-2.43239e-08,1.65284e-11,28393.3,12.5062], Tmin=(100,'K'), Tmax=(929.088,'K')), NASAPolynomial(coeffs=[18.6875,0.00574194,-9.61607e-08,-4.92876e-11,1.62459e-16,23761.4,-77.9814], Tmin=(929.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
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
    label = 'OH(D)(132)',
    structure = SMILES('[OH]'),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92814e-05,-5.32177e-07,1.01951e-09,-3.85951e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604011,-1.39759e-08,-2.13452e-11,2.4807e-15,3579.39,4.57799], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CCC=C=C[O](22612)',
    structure = SMILES('C#CCC=C=C[O]'),
    E0 = (331.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.36274,'amu*angstrom^2'), symmetry=1, barrier=(31.332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36681,'amu*angstrom^2'), symmetry=1, barrier=(31.4257,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02349,0.0552307,-3.16793e-05,-6.99827e-09,9.08567e-12,39972.1,22.8024], Tmin=(100,'K'), Tmax=(951.932,'K')), NASAPolynomial(coeffs=[16.9064,0.0130554,-3.9288e-06,6.74845e-10,-4.85801e-14,35835.2,-58.8869], Tmin=(951.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=C[CH]C=CO(25690)',
    structure = SMILES('[O]C=C=C[CH]C=CO'),
    E0 = (54.9308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277228,0.0555203,2.61182e-05,-1.00683e-07,5.07216e-11,6765.37,27.8461], Tmin=(100,'K'), Tmax=(908.931,'K')), NASAPolynomial(coeffs=[29.3612,-0.0024519,6.23615e-06,-1.34677e-09,8.77881e-14,-1414.07,-125.601], Tmin=(908.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.9308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CCC=[C]O(25687)',
    structure = SMILES('[O]C=C=CCC=[C]O'),
    E0 = (192.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0116,'amu*angstrom^2'), symmetry=1, barrier=(23.2586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00756,'amu*angstrom^2'), symmetry=1, barrier=(23.1659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00711,'amu*angstrom^2'), symmetry=1, barrier=(23.1554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169504,0.0696731,-3.85264e-05,-1.67955e-08,1.63885e-11,23357.9,31.6896], Tmin=(100,'K'), Tmax=(926.396,'K')), NASAPolynomial(coeffs=[22.9819,0.008834,-9.96421e-07,7.94876e-11,-7.81948e-15,17515.2,-85.3453], Tmin=(926.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'OC#C[CH]C[C]=CO(29143)',
    structure = SMILES('OC#C[CH]C[C]=CO'),
    E0 = (263.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2100,2250,500,550,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.428136,0.0834538,-9.25546e-05,5.03147e-08,-1.02696e-11,31874.8,32.2805], Tmin=(100,'K'), Tmax=(1361.95,'K')), NASAPolynomial(coeffs=[20.6849,0.0109734,-1.19384e-06,-5.15847e-11,1.19414e-14,27095.1,-72.5344], Tmin=(1361.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=[C]CC=CO(25692)',
    structure = SMILES('[O]C=C=[C]CC=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=CCC=C[O](22588)',
    structure = SMILES('[O]C=C=CCC=C[O]'),
    E0 = (94.6684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24294,'amu*angstrom^2'), symmetry=1, barrier=(28.5775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24272,'amu*angstrom^2'), symmetry=1, barrier=(28.5727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4556.79,'J/mol'), sigma=(7.12408,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=711.76 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371633,0.0609434,-8.25709e-06,-4.93002e-08,2.78157e-11,11534,29.2454], Tmin=(100,'K'), Tmax=(935.092,'K')), NASAPolynomial(coeffs=[23.8145,0.00848651,-8.25121e-07,9.47677e-11,-1.27022e-14,5058.87,-93.4597], Tmin=(935.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'O[CH]C#CC[C]=CO(29144)',
    structure = SMILES('O[CH]C#CC[C]=CO'),
    E0 = (272.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,1685,370,2100,2250,500,550,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0626319,0.0816345,-8.61059e-05,3.72226e-08,-2.70289e-12,32879.9,31.1856], Tmin=(100,'K'), Tmax=(859.877,'K')), NASAPolynomial(coeffs=[19.141,0.0133617,-2.72823e-06,2.73208e-10,-1.19288e-14,28841.9,-62.378], Tmin=(859.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'O=C=[C][CH]CC=CO(25695)',
    structure = SMILES('O=C=[C][CH]CC=CO'),
    E0 = (119.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05798,'amu*angstrom^2'), symmetry=1, barrier=(24.3251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05848,'amu*angstrom^2'), symmetry=1, barrier=(24.3365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.058,'amu*angstrom^2'), symmetry=1, barrier=(24.3256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0596,'amu*angstrom^2'), symmetry=1, barrier=(24.3624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377521,0.0817048,-8.38588e-05,4.18731e-08,-7.97521e-12,14521.9,31.8005], Tmin=(100,'K'), Tmax=(1386.82,'K')), NASAPolynomial(coeffs=[21.854,0.0124194,-3.33443e-06,4.7915e-10,-2.92203e-14,8852.24,-80.9334], Tmin=(1386.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'OC=[C]C=C[C]=CO(29145)',
    structure = SMILES('OC=[C]C=C[C]=CO'),
    E0 = (126.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.84374,0.109539,-0.000129866,6.8666e-08,-1.28873e-11,15448,33.4179], Tmin=(100,'K'), Tmax=(1611.54,'K')), NASAPolynomial(coeffs=[29.4599,-0.00516106,9.02421e-06,-2.08208e-09,1.49869e-13,9518.61,-123.937], Tmin=(1611.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'O[C]=[C]CC=C=CO(29146)',
    structure = SMILES('O[C]=[C]CC=C=CO'),
    E0 = (289.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.956189,'amu*angstrom^2'), symmetry=1, barrier=(21.9847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952474,'amu*angstrom^2'), symmetry=1, barrier=(21.8992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957255,'amu*angstrom^2'), symmetry=1, barrier=(22.0092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955655,'amu*angstrom^2'), symmetry=1, barrier=(21.9724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17587,0.0918469,-0.000103785,5.44468e-08,-1.04981e-11,35003.5,35.5826], Tmin=(100,'K'), Tmax=(1480.58,'K')), NASAPolynomial(coeffs=[25.5528,0.00345682,2.1557e-06,-6.35668e-10,4.89189e-14,28862,-97.8687], Tmin=(1480.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]CC=C=CO(25694)',
    structure = SMILES('[O]C=[C]CC=C=CO'),
    E0 = (191.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11134,'amu*angstrom^2'), symmetry=1, barrier=(25.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11128,'amu*angstrom^2'), symmetry=1, barrier=(25.5506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11143,'amu*angstrom^2'), symmetry=1, barrier=(25.554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0688376,0.0724563,-3.61978e-05,-2.66638e-08,2.18221e-11,23140.1,29.8906], Tmin=(100,'K'), Tmax=(916.915,'K')), NASAPolynomial(coeffs=[25.9182,0.00423282,1.55867e-06,-4.19886e-10,2.63005e-14,16476.8,-103.576], Tmin=(916.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]CC=C=C[O](22626)',
    structure = SMILES('[CH]=[C]CC=C=C[O]'),
    E0 = (646.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06585,'amu*angstrom^2'), symmetry=1, barrier=(24.506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06786,'amu*angstrom^2'), symmetry=1, barrier=(24.5523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3942.94,'J/mol'), sigma=(6.38069,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.88 K, Pc=34.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874198,0.0604722,-4.9312e-05,1.32472e-08,1.2772e-12,77928.4,27.0485], Tmin=(100,'K'), Tmax=(985.722,'K')), NASAPolynomial(coeffs=[16.5106,0.0141667,-4.93912e-06,8.83067e-10,-6.24589e-14,74012.8,-52.3884], Tmin=(985.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C=[C]CC=C=C[O](25698)',
    structure = SMILES('[O]C=[C]CC=C=C[O]'),
    E0 = (332.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2191,'amu*angstrom^2'), symmetry=1, barrier=(28.0296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22259,'amu*angstrom^2'), symmetry=1, barrier=(28.1097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342925,0.0654961,-3.29055e-05,-2.01592e-08,1.68742e-11,40137.1,29.822], Tmin=(100,'K'), Tmax=(936.243,'K')), NASAPolynomial(coeffs=[22.606,0.00801836,-1.12086e-06,1.48001e-10,-1.43778e-14,34318.8,-84.9277], Tmin=(936.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=CO(18753)',
    structure = SMILES('[CH2][C]=CO'),
    E0 = (186.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23523,'amu*angstrom^2'), symmetry=1, barrier=(28.4004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2351,'amu*angstrom^2'), symmetry=1, barrier=(28.3973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24497,0.0260528,1.3484e-05,-5.00525e-08,2.54383e-11,22526.5,14.0801], Tmin=(100,'K'), Tmax=(898.827,'K')), NASAPolynomial(coeffs=[16.2027,-0.00210248,3.79693e-06,-8.3211e-10,5.63273e-14,18645.6,-59.4], Tmin=(898.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O]C=[C]C=C[C]=CO(29147)',
    structure = SMILES('[O]C=[C]C=C[C]=CO'),
    E0 = (267.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48247,'amu*angstrom^2'), symmetry=1, barrier=(34.0848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4964,'amu*angstrom^2'), symmetry=1, barrier=(34.4052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49353,'amu*angstrom^2'), symmetry=1, barrier=(34.3392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17784,0.0995255,-0.000115587,6.03866e-08,-1.12584e-11,32434,33.134], Tmin=(100,'K'), Tmax=(1612.68,'K')), NASAPolynomial(coeffs=[27.1561,-0.00276498,7.02459e-06,-1.65532e-09,1.19747e-13,26813,-110.497], Tmin=(1612.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=C=[C]C[C]=CO(29148)',
    structure = SMILES('[O]C=C=[C]C[C]=CO'),
    E0 = (428.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10287,'amu*angstrom^2'), symmetry=1, barrier=(25.3571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10245,'amu*angstrom^2'), symmetry=1, barrier=(25.3474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10246,'amu*angstrom^2'), symmetry=1, barrier=(25.3477,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01625,0.0877556,-9.81082e-05,5.05072e-08,-9.54698e-12,51783.5,33.7698], Tmin=(100,'K'), Tmax=(1510.84,'K')), NASAPolynomial(coeffs=[25.2863,0.00244976,2.14207e-06,-5.93198e-10,4.45036e-14,45624,-98.0607], Tmin=(1510.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=CC[C]=[C]O(29149)',
    structure = SMILES('[O]C=C=CC[C]=[C]O'),
    E0 = (430.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.989153,'amu*angstrom^2'), symmetry=1, barrier=(22.7426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990709,'amu*angstrom^2'), symmetry=1, barrier=(22.7783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989564,'amu*angstrom^2'), symmetry=1, barrier=(22.752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.510351,0.0818406,-8.95382e-05,4.62166e-08,-8.89101e-12,51989.5,34.6068], Tmin=(100,'K'), Tmax=(1444.96,'K')), NASAPolynomial(coeffs=[22.9358,0.00629248,-6.32667e-08,-1.6251e-10,1.52652e-14,46324.9,-83.2923], Tmin=(1444.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'O=C=[C][CH]C[C]=CO(29150)',
    structure = SMILES('O=C=[C][CH]C[C]=CO'),
    E0 = (357.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1670,1700,300,440,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00937,'amu*angstrom^2'), symmetry=1, barrier=(23.2073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00958,'amu*angstrom^2'), symmetry=1, barrier=(23.2122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00955,'amu*angstrom^2'), symmetry=1, barrier=(23.2115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00917,'amu*angstrom^2'), symmetry=1, barrier=(23.2027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.030048,0.0819375,-9.39217e-05,5.28016e-08,-1.14314e-11,43108.6,31.0197], Tmin=(100,'K'), Tmax=(1142.36,'K')), NASAPolynomial(coeffs=[19.3464,0.0140907,-4.83436e-06,8.11626e-10,-5.37277e-14,38681.6,-65.0395], Tmin=(1142.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCJC(C)=C=O) + radical(Cds_S) + radical(CCCJ=C=O)"""),
)

species(
    label = 'OC=[C]CC=C1[CH]O1(29151)',
    structure = SMILES('OC=[C]CC=C1[CH]O1'),
    E0 = (217.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311273,0.0517352,3.94608e-05,-1.17245e-07,5.70562e-11,26287,27.108], Tmin=(100,'K'), Tmax=(913.94,'K')), NASAPolynomial(coeffs=[31.2318,-0.00628158,7.79276e-06,-1.58732e-09,1.00963e-13,17406.2,-136.948], Tmin=(913.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C1[CH]CC1=CO(29132)',
    structure = SMILES('[O]C=C1[CH]CC1=CO'),
    E0 = (51.1483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580919,0.0459017,5.30866e-05,-1.26455e-07,5.89378e-11,6302.11,22.9576], Tmin=(100,'K'), Tmax=(917.612,'K')), NASAPolynomial(coeffs=[28.9824,-0.00147031,5.57963e-06,-1.1648e-09,7.14797e-14,-2128.11,-129.155], Tmin=(917.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.1483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(12methylenecyclobutane) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = 'OC=[C]CC1[C]=CO1(29152)',
    structure = SMILES('OC=[C]CC1[C]=CO1'),
    E0 = (309.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0890768,0.0572281,2.70846e-05,-1.06177e-07,5.34324e-11,37405.9,27.5852], Tmin=(100,'K'), Tmax=(914.748,'K')), NASAPolynomial(coeffs=[32.052,-0.00677526,7.80019e-06,-1.57854e-09,1.00307e-13,28388.5,-141.095], Tmin=(914.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]=CCC1=CO(29153)',
    structure = SMILES('[O]C1[C]=CCC1=CO'),
    E0 = (216.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991361,0.0459413,2.37909e-05,-7.43341e-08,3.44641e-11,26179.2,27.3261], Tmin=(100,'K'), Tmax=(948.875,'K')), NASAPolynomial(coeffs=[20.9451,0.0123389,-2.94157e-06,5.50644e-10,-4.7185e-14,20118.5,-79.8747], Tmin=(948.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(4-Methylenecyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'OC=C=CC=C=CO(29154)',
    structure = SMILES('OC=C=CC=C=CO'),
    E0 = (-42.5067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.79511,0.104175,-0.000117247,5.85778e-08,-1.04185e-11,-4828.48,32.4729], Tmin=(100,'K'), Tmax=(1696.98,'K')), NASAPolynomial(coeffs=[29.3807,-0.00417545,7.2608e-06,-1.62395e-09,1.13479e-13,-11068.1,-125.982], Tmin=(1696.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.5067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=C=C=CCC=CO(25707)',
    structure = SMILES('O=C=C=CCC=CO'),
    E0 = (-2.69765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23642,0.0443849,6.55999e-06,-5.33988e-08,2.74403e-11,-209.807,11.9333], Tmin=(100,'K'), Tmax=(929.65,'K')), NASAPolynomial(coeffs=[19.8995,0.00620394,2.03163e-07,-1.03377e-10,1.90958e-15,-5499.98,-86.5333], Tmin=(929.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.69765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O]C=[C]C[CH][C]=CO(29155)',
    structure = SMILES('[O]C=[C]C[CH][C]=CO'),
    E0 = (404.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3025,407.5,1350,352.5,218.578,219.288,219.358,220.606],'cm^-1')),
        HinderedRotor(inertia=(0.662305,'amu*angstrom^2'), symmetry=1, barrier=(22.4896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.67626,'amu*angstrom^2'), symmetry=1, barrier=(22.5131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.676313,'amu*angstrom^2'), symmetry=1, barrier=(22.5182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668726,'amu*angstrom^2'), symmetry=1, barrier=(22.5044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.000353015,0.0722127,-3.94057e-05,-2.02052e-08,1.8683e-11,48793.6,30.3692], Tmin=(100,'K'), Tmax=(921.012,'K')), NASAPolynomial(coeffs=[24.6571,0.00648526,2.78239e-07,-1.70267e-10,9.31855e-15,42497.4,-96.0835], Tmin=(921.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(Cds_S) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C[C]C[C]=CO(29156)',
    structure = SMILES('[O]C=C[C]C[C]=CO'),
    E0 = (473.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60944,0.0919313,-9.76574e-05,4.76631e-08,-8.48242e-12,57194.1,36.708], Tmin=(100,'K'), Tmax=(1639.61,'K')), NASAPolynomial(coeffs=[26.4523,0.00298311,2.46103e-06,-6.66545e-10,4.87396e-14,50746,-104.151], Tmin=(1639.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=C[CH]C[C]=CO(29157)',
    structure = SMILES('[O][C]=C[CH]C[C]=CO'),
    E0 = (406.276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3025,407.5,1350,352.5,369.691,369.691,369.691,369.691],'cm^-1')),
        HinderedRotor(inertia=(0.201904,'amu*angstrom^2'), symmetry=1, barrier=(19.5817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201904,'amu*angstrom^2'), symmetry=1, barrier=(19.5817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201904,'amu*angstrom^2'), symmetry=1, barrier=(19.5817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201904,'amu*angstrom^2'), symmetry=1, barrier=(19.5817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.236271,0.0694491,-4.1798e-05,-1.0265e-08,1.32253e-11,49011.5,32.1744], Tmin=(100,'K'), Tmax=(932.785,'K')), NASAPolynomial(coeffs=[21.7328,0.0110662,-2.26516e-06,3.26346e-10,-2.45725e-14,43530.8,-77.9199], Tmin=(932.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C[CH]C=[C][CH]O(29158)',
    structure = SMILES('[O]C=C[CH]C=[C][CH]O'),
    E0 = (289.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,446.538,447.288,447.375,447.396],'cm^-1')),
        HinderedRotor(inertia=(0.000843895,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240437,'amu*angstrom^2'), symmetry=1, barrier=(34.0861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240208,'amu*angstrom^2'), symmetry=1, barrier=(34.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240007,'amu*angstrom^2'), symmetry=1, barrier=(34.0838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826411,0.0491241,2.07576e-05,-7.64109e-08,3.69046e-11,34906.5,31.7859], Tmin=(100,'K'), Tmax=(931.012,'K')), NASAPolynomial(coeffs=[22.2974,0.00973009,-9.28196e-07,9.45643e-11,-1.28841e-14,28617.9,-82.5663], Tmin=(931.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJO) + radical(C=CCJC=C) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C[C]=C[CH][C]=CO(29159)',
    structure = SMILES('[O]C[C]=C[CH][C]=CO'),
    E0 = (493.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,408.169,408.205,408.228,408.237],'cm^-1')),
        HinderedRotor(inertia=(0.001012,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138778,'amu*angstrom^2'), symmetry=1, barrier=(16.4049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138766,'amu*angstrom^2'), symmetry=1, barrier=(16.4054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649544,'amu*angstrom^2'), symmetry=1, barrier=(76.8158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5367,0.0690126,-5.9802e-05,2.24882e-08,-1.82206e-12,59534.6,32.1516], Tmin=(100,'K'), Tmax=(984.07,'K')), NASAPolynomial(coeffs=[15.9092,0.0204907,-7.12539e-06,1.22138e-09,-8.24717e-14,55833,-45.2004], Tmin=(984.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCOJ) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C(C=C=C[O])=CO(25407)',
    structure = SMILES('[CH2]C(C=C=C[O])=CO'),
    E0 = (72.2858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35624,'amu*angstrom^2'), symmetry=1, barrier=(31.1827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35316,'amu*angstrom^2'), symmetry=1, barrier=(31.1118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35546,'amu*angstrom^2'), symmetry=1, barrier=(31.1646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42989,0.0741445,-1.93668e-05,-5.92742e-08,3.73878e-11,8875.35,26.7673], Tmin=(100,'K'), Tmax=(904.884,'K')), NASAPolynomial(coeffs=[31.7552,-0.00493166,6.95647e-06,-1.48733e-09,9.9297e-14,463.293,-139.586], Tmin=(904.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.2858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=C[O](22438)',
    structure = SMILES('C=C[C]=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61747,'amu*angstrom^2'), symmetry=1, barrier=(37.1889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98077,0.0322088,7.53766e-06,-4.54875e-08,2.38742e-11,27216,16.1107], Tmin=(100,'K'), Tmax=(899.949,'K')), NASAPolynomial(coeffs=[16.1069,0.00249504,1.93925e-06,-5.05324e-10,3.47549e-14,23334.1,-57.9915], Tmin=(899.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[C]=CO(27807)',
    structure = SMILES('[C]=CO'),
    E0 = (391.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.23325,'amu*angstrom^2'), symmetry=1, barrier=(28.3547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88904,0.0147275,1.32236e-05,-4.12494e-08,2.11475e-11,47130.8,9.58665], Tmin=(100,'K'), Tmax=(884.362,'K')), NASAPolynomial(coeffs=[13.839,-0.00784764,5.80008e-06,-1.19218e-09,8.19757e-14,44140.2,-47.8537], Tmin=(884.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
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
    label = 'C#C[CH]C[C]=CO(27789)',
    structure = SMILES('C#C[CH]C[C]=CO'),
    E0 = (405.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2175,525,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.997839,'amu*angstrom^2'), symmetry=1, barrier=(22.9423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.997839,'amu*angstrom^2'), symmetry=1, barrier=(22.9423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.55921,'amu*angstrom^2'), symmetry=1, barrier=(58.8413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.997112,'amu*angstrom^2'), symmetry=1, barrier=(22.9256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.259851,0.0766703,-8.43788e-05,4.52067e-08,-8.94628e-12,48885.8,28.3377], Tmin=(100,'K'), Tmax=(1459.4,'K')), NASAPolynomial(coeffs=[19.3254,0.00848098,6.19949e-07,-4.33581e-10,3.8988e-14,44714.3,-68.2601], Tmin=(1459.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'O=CC#CC[C]=CO(29160)',
    structure = SMILES('O=CC#CC[C]=CO'),
    E0 = (157.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1685,370,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.988434,'amu*angstrom^2'), symmetry=1, barrier=(22.726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.987841,'amu*angstrom^2'), symmetry=1, barrier=(22.7124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.987813,'amu*angstrom^2'), symmetry=1, barrier=(22.7118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.987422,'amu*angstrom^2'), symmetry=1, barrier=(22.7028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0767663,0.0742406,-7.69344e-05,3.89083e-08,-7.53118e-12,19060.5,29.2364], Tmin=(100,'K'), Tmax=(1342,'K')), NASAPolynomial(coeffs=[19.9123,0.0117738,-3.37494e-06,5.09071e-10,-3.18762e-14,14037.8,-71.1714], Tmin=(1342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S)"""),
)

species(
    label = 'O=CC=[C]C[C]=CO(29161)',
    structure = SMILES('O=CC=[C]C[C]=CO'),
    E0 = (232.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.782466,'amu*angstrom^2'), symmetry=1, barrier=(17.9904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782309,'amu*angstrom^2'), symmetry=1, barrier=(17.9868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782422,'amu*angstrom^2'), symmetry=1, barrier=(17.9894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.780185,'amu*angstrom^2'), symmetry=1, barrier=(17.938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272915,0.0784276,-8.23303e-05,4.2882e-08,-8.78412e-12,28058.7,30.087], Tmin=(100,'K'), Tmax=(1187.74,'K')), NASAPolynomial(coeffs=[17.5097,0.0203784,-9.01934e-06,1.73298e-09,-1.2287e-13,23964.2,-56.0359], Tmin=(1187.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'O=C=C[CH]C[C]=CO(29162)',
    structure = SMILES('O=C=C[CH]C[C]=CO'),
    E0 = (154.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0439598,0.0812126,-8.76673e-05,4.70442e-08,-9.78041e-12,18781.1,31.3614], Tmin=(100,'K'), Tmax=(1185.37,'K')), NASAPolynomial(coeffs=[19.1367,0.0164879,-5.76275e-06,9.79968e-10,-6.52357e-14,14233.9,-64.4357], Tmin=(1185.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[O]C=C[CH]C=C=CO(25693)',
    structure = SMILES('[O]C=C[CH]C=C=CO'),
    E0 = (54.9308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277228,0.0555203,2.61182e-05,-1.00683e-07,5.07216e-11,6765.37,27.8461], Tmin=(100,'K'), Tmax=(908.931,'K')), NASAPolynomial(coeffs=[29.3612,-0.0024519,6.23615e-06,-1.34677e-09,8.77881e-14,-1414.07,-125.601], Tmin=(908.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.9308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = 'O=CC=CC[C]=[C]O(29163)',
    structure = SMILES('O=CC=CC[C]=[C]O'),
    E0 = (234.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.635524,'amu*angstrom^2'), symmetry=1, barrier=(14.612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635343,'amu*angstrom^2'), symmetry=1, barrier=(14.6078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635498,'amu*angstrom^2'), symmetry=1, barrier=(14.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635505,'amu*angstrom^2'), symmetry=1, barrier=(14.6115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739949,0.0728951,-7.48177e-05,3.96877e-08,-8.51479e-12,28266.8,31.0694], Tmin=(100,'K'), Tmax=(1116.68,'K')), NASAPolynomial(coeffs=[13.7103,0.0264344,-1.24081e-05,2.42836e-09,-1.73189e-13,25370,-32.9362], Tmin=(1116.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]CC=CC=O(25734)',
    structure = SMILES('[O]C=[C]CC=CC=O'),
    E0 = (135.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.795126,'amu*angstrom^2'), symmetry=1, barrier=(18.2815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795313,'amu*angstrom^2'), symmetry=1, barrier=(18.2858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794299,'amu*angstrom^2'), symmetry=1, barrier=(18.2625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582611,0.068327,-5.86277e-05,2.45842e-08,-4.08711e-12,16458.3,29.9183], Tmin=(100,'K'), Tmax=(1436.45,'K')), NASAPolynomial(coeffs=[17.4345,0.0214009,-9.62585e-06,1.84222e-09,-1.29125e-13,11616.9,-57.4857], Tmin=(1436.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'OC=C1CC=[C][CH]O1(29164)',
    structure = SMILES('OC=C1CC=[C][CH]O1'),
    E0 = (63.5399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1649,0.0334323,7.73013e-05,-1.35077e-07,5.63661e-11,7770.47,21.6882], Tmin=(100,'K'), Tmax=(954.606,'K')), NASAPolynomial(coeffs=[23.6819,0.0117555,-2.83225e-06,6.36104e-10,-6.11456e-14,160.215,-103.241], Tmin=(954.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.5399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'O=CC=CC=C=CO(25728)',
    structure = SMILES('O=CC=CC=C=CO'),
    E0 = (-97.7808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0434416,0.0759566,-5.83246e-05,7.98723e-09,5.69611e-12,-11603.1,26.3036], Tmin=(100,'K'), Tmax=(979.023,'K')), NASAPolynomial(coeffs=[22.5942,0.0119884,-4.01707e-06,7.64445e-10,-5.83952e-14,-17402.6,-89.4112], Tmin=(979.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.7808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=CC1=CCC1=CO(29123)',
    structure = SMILES('O=CC1=CCC1=CO'),
    E0 = (-103.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32454,0.0544497,2.4609e-05,-9.34728e-08,4.57186e-11,-12234,22.8026], Tmin=(100,'K'), Tmax=(932.608,'K')), NASAPolynomial(coeffs=[28.6901,0.000263098,3.23675e-06,-6.1642e-10,3.08926e-14,-20459.1,-127.797], Tmin=(932.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    label = '[C]=CC[C]=CO(28208)',
    structure = SMILES('[C]=CC[C]=CO'),
    E0 = (675.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.834503,'amu*angstrom^2'), symmetry=1, barrier=(19.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834445,'amu*angstrom^2'), symmetry=1, barrier=(19.1855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835194,'amu*angstrom^2'), symmetry=1, barrier=(19.2028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597635,0.0609935,-6.39407e-05,3.20564e-08,-6.01275e-12,81426.4,25.7518], Tmin=(100,'K'), Tmax=(1483.52,'K')), NASAPolynomial(coeffs=[17.8834,0.00639117,-6.48154e-07,-1.88037e-11,4.66406e-15,77177.5,-61.4945], Tmin=(1483.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(675.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
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
    E0 = (191.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (274.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (321.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (416.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (449.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (260.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (359.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (347.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (425.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (441.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (349.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (341.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (587.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (283.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (282.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (466.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (308.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (675.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (544.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (456.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (479.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (640.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (642.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (568.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (378.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (317.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (316.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (247.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (208.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (230.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (432.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (495.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (429.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (352.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (518.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (361.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (651.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (811.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (378.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (287.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (428.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (353.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (367.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (430.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (381.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (249.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (284.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (199.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (743.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['C=C=CO(12571)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=[C]C1CC1=CO(29139)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[O]C=C=CC=C=CO(29140)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.47e+08,'cm^3/(mol*s)'), n=1.64, Ea=(10.711,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2592 used for Cds-CdH_Ca;HJ
Exact match found for rate rule [Cds-CdH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]C=C=CCC#CO(29141)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'O=C=C=CC[C]=CO(29142)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=C[O](8556)', 'C=C=CO(12571)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Ca;CJ] for rate rule [Cds-HH_Ca;CdsJ=Cdd]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', 'C#CCC=C=C[O](22612)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.508e+07,'cm^3/(mol*s)'), n=1.628, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 211 used for Ct-H_Ct-Cs;OJ_pri
Exact match found for rate rule [Ct-H_Ct-Cs;OJ_pri]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C=C[CH]C=CO(25690)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C=CCC=[C]O(25687)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC#C[CH]C[C]=CO(29143)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C=[C]CC=CO(25692)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.85004e+08,'s^-1'), n=1.20667, Ea=(158.922,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] + [R3H_SS_Cs;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C=CCC=C[O](22588)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O[CH]C#CC[C]=CO(29144)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['O=C=[C][CH]CC=CO(25695)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.08955e+06,'s^-1'), n=1.58778, Ea=(92.3497,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R5H;Cd_rad_out;XH_out] for rate rule [R5H;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['OC=[C]C=C[C]=CO(29145)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(492144,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O[C]=[C]CC=C=CO(29146)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;XH_out] for rate rule [R7HJ_1;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=[C]CC=C=CO(25694)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(40268.4,'s^-1'), n=2.03024, Ea=(116.997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [RnH;O_rad_out;XH_out] + [R8Hall;Y_rad_out;XH_out] for rate rule [R8Hall;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OH(D)(132)', '[CH]=[C]CC=C=C[O](22626)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[O]C=[C]CC=C=C[O](25698)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=C[O](8556)', '[CH2][C]=CO(18753)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[O]C=[C]C=C[C]=CO(29147)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[O]C=C=[C]C[C]=CO(29148)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[O]C=C=CC[C]=[C]O(29149)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', 'O=C=[C][CH]C[C]=CO(29150)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['OC=[C]CC=C1[CH]O1(29151)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C1[CH]CC1=CO(29132)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['OC=[C]CC1[C]=CO1(29152)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C1[C]=CCC1=CO(29153)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.01659e+11,'s^-1'), n=0.0800815, Ea=(56.4798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra;radadd_intra_cddouble] + [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['OC=C=CC=C=CO(29154)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['O=C=C=CCC=CO(25707)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C=[C]C[CH][C]=CO(29155)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C=C[C]C[C]=CO(29156)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O][C]=C[CH]C[C]=CO(29157)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C=C[CH]C=[C][CH]O(29158)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C[C]=C[CH][C]=CO(29159)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[C]=C[O](22438)', '[C]=CO(27807)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(T)(63)', 'C#C[CH]C[C]=CO(27789)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(8)', 'O=CC#CC[C]=CO(29160)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]=CO(18753)', 'C#CC=O(21959)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O=CC=[C]C[C]=CO(29161)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['O=C=C[CH]C[C]=CO(29162)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=C[CH]C=C=CO(25693)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.65652e+07,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O=CC=CC[C]=[C]O(29163)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleDe] for rate rule [R5HJ_1;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['[O]C=[C]CC=CC=O(25734)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.48104e+08,'s^-1'), n=1.45648, Ea=(190.329,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe;XH_out] for rate rule [R6HJ_3;Cd_rad_out_singleDe;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['OC=C1CC=[C][CH]O1(29164)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.19745e+11,'s^-1'), n=0.440371, Ea=(58.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['O=CC=CC=C=CO(25728)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]C=C=CC[C]=CO(25689)'],
    products = ['O=CC1=CCC1=CO(29123)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=O(373)', '[C]=CC[C]=CO(28208)'],
    products = ['[O]C=C=CC[C]=CO(25689)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '5261',
    isomers = [
        '[O]C=C=CC[C]=CO(25689)',
    ],
    reactants = [
        ('C=C=CO(12571)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5261',
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


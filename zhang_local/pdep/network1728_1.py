species(
    label = '[CH2]C([CH2])C[CH]CC(261)',
    structure = SMILES('[CH2]C([CH2])C[CH]CC'),
    E0 = (383.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,357.527,506.612,3176.78],'cm^-1')),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0275299,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791632,0.0726375,-5.49107e-05,2.77146e-08,-6.51272e-12,46212.1,33.0159], Tmin=(100,'K'), Tmax=(961.311,'K')), NASAPolynomial(coeffs=[6.14615,0.0503581,-2.01475e-05,3.60712e-09,-2.43475e-13,45182.6,7.39456], Tmin=(961.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC(36)',
    structure = SMILES('C=CCC'),
    E0 = (-16.5218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,385.42],'cm^-1')),
        HinderedRotor(inertia=(0.10941,'amu*angstrom^2'), symmetry=1, barrier=(11.0707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106381,'amu*angstrom^2'), symmetry=1, barrier=(11.1092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64174,0.0203702,3.0581e-05,-4.85559e-08,1.85232e-11,-1929.81,14.9537], Tmin=(100,'K'), Tmax=(991.794,'K')), NASAPolynomial(coeffs=[7.80835,0.0229241,-8.65887e-06,1.6005e-09,-1.13877e-13,-4105.1,-15.7294], Tmin=(991.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.5218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C=C(87)',
    structure = SMILES('[CH2]C=C'),
    E0 = (157.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.570287,'amu*angstrom^2'), symmetry=1, barrier=(32.8573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3193,0.00566487,4.27449e-05,-5.78831e-08,2.21699e-11,18990.6,9.19646], Tmin=(100,'K'), Tmax=(951.999,'K')), NASAPolynomial(coeffs=[7.55715,0.0114811,-3.63952e-06,6.63584e-10,-4.95318e-14,17113.3,-16.6624], Tmin=(951.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=C)C[CH]CC(6234)',
    structure = SMILES('[CH2]C(=C)C[CH]CC'),
    E0 = (242.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839832,0.0669513,-3.9312e-05,1.16119e-08,-1.42369e-12,29319.6,28.8165], Tmin=(100,'K'), Tmax=(1779.23,'K')), NASAPolynomial(coeffs=[13.2578,0.0390338,-1.5776e-05,2.79316e-09,-1.84568e-13,24900.7,-38.2481], Tmin=(1779.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C([CH2])C=CCC(6235)',
    structure = SMILES('[CH2]C([CH2])C=CCC'),
    E0 = (303.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,295.741,2046.92],'cm^-1')),
        HinderedRotor(inertia=(0.0019129,'amu*angstrom^2'), symmetry=1, barrier=(0.11963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167398,'amu*angstrom^2'), symmetry=1, barrier=(10.3981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166943,'amu*angstrom^2'), symmetry=1, barrier=(10.399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00191243,'amu*angstrom^2'), symmetry=1, barrier=(0.119635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1736,'amu*angstrom^2'), symmetry=1, barrier=(72.9413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665302,0.0641745,-2.33028e-05,-1.24013e-08,9.1857e-12,36597.1,31.6344], Tmin=(100,'K'), Tmax=(972.579,'K')), NASAPolynomial(coeffs=[12.257,0.037008,-1.3033e-05,2.23951e-09,-1.50614e-13,33372.4,-28.953], Tmin=(972.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])CC=CC(6236)',
    structure = SMILES('[CH2]C([CH2])CC=CC'),
    E0 = (301.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,357.979,1864.27],'cm^-1')),
        HinderedRotor(inertia=(0.0702683,'amu*angstrom^2'), symmetry=1, barrier=(6.39097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702646,'amu*angstrom^2'), symmetry=1, barrier=(6.39098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702652,'amu*angstrom^2'), symmetry=1, barrier=(6.39147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702761,'amu*angstrom^2'), symmetry=1, barrier=(6.39086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75564,'amu*angstrom^2'), symmetry=1, barrier=(68.726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341733,0.0710699,-4.85002e-05,1.8387e-08,-2.91057e-12,36343.3,32.1028], Tmin=(100,'K'), Tmax=(1467.29,'K')), NASAPolynomial(coeffs=[13.2106,0.035987,-1.26341e-05,2.09073e-09,-1.33892e-13,32566.9,-34.9158], Tmin=(1467.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = 'C=CC[CH]CC(216)',
    structure = SMILES('C=CC[CH]CC'),
    E0 = (130.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,243.588,1471.94,1472.19],'cm^-1')),
        HinderedRotor(inertia=(0.00283977,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18429,'amu*angstrom^2'), symmetry=1, barrier=(7.7782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18492,'amu*angstrom^2'), symmetry=1, barrier=(7.7775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184441,'amu*angstrom^2'), symmetry=1, barrier=(7.77788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47064,0.051542,-2.30987e-05,3.77216e-09,-2.75739e-14,15775.1,26.2413], Tmin=(100,'K'), Tmax=(1745.76,'K')), NASAPolynomial(coeffs=[12.9328,0.0324594,-1.2872e-05,2.22284e-09,-1.43099e-13,10678.9,-38.5776], Tmin=(1745.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CC(39)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1553.58],'cm^-1')),
        HinderedRotor(inertia=(0.00260974,'amu*angstrom^2'), symmetry=1, barrier=(4.49727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19501,'amu*angstrom^2'), symmetry=1, barrier=(4.48366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194146,'amu*angstrom^2'), symmetry=1, barrier=(4.4638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59221,0.0286125,-6.18652e-06,-3.09336e-09,1.2171e-12,30768.7,19.1946], Tmin=(100,'K'), Tmax=(1492.38,'K')), NASAPolynomial(coeffs=[6.22741,0.0257398,-1.0205e-05,1.78668e-09,-1.17174e-13,28918.5,-2.36199], Tmin=(1492.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH2](6136)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00132719,'amu*angstrom^2'), symmetry=1, barrier=(2.41051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132749,'amu*angstrom^2'), symmetry=1, barrier=(2.41088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34321,0.013802,2.16426e-06,-5.76329e-09,1.61332e-12,58271,14.955], Tmin=(100,'K'), Tmax=(1447.11,'K')), NASAPolynomial(coeffs=[4.39505,0.0167645,-6.99091e-06,1.25741e-09,-8.38108e-14,57351.9,7.36811], Tmin=(1447.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH3](11)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1408.13,1408.49,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48742e-06,-3.32748e-09,8.49957e-13,16285.6,0.351741], Tmin=(100,'K'), Tmax=(1337.63,'K')), NASAPolynomial(coeffs=[3.54146,0.00476787,-1.82148e-06,3.28877e-10,-2.22546e-14,16224,1.66035], Tmin=(1337.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH2]C([CH2])CC=C(6127)',
    structure = SMILES('[CH2]C([CH2])CC=C'),
    E0 = (337.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,284.918,1590.2],'cm^-1')),
        HinderedRotor(inertia=(0.00207914,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00207568,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135159,'amu*angstrom^2'), symmetry=1, barrier=(7.78346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21238,'amu*angstrom^2'), symmetry=1, barrier=(69.776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29623,0.0520097,-1.6468e-05,-1.38423e-08,9.25212e-12,40640.7,26.8619], Tmin=(100,'K'), Tmax=(947.416,'K')), NASAPolynomial(coeffs=[10.6915,0.0302258,-1.0292e-05,1.73515e-09,-1.15664e-13,38057.9,-22.1931], Tmin=(947.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C[CH]CC(258)',
    structure = SMILES('[CH2][C](C)C[CH]CC'),
    E0 = (363.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,273.042,1296.09,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126332,'amu*angstrom^2'), symmetry=1, barrier=(0.492855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773309,0.0795277,-9.41071e-05,8.49346e-08,-3.20446e-11,43841.8,33.1201], Tmin=(100,'K'), Tmax=(837.895,'K')), NASAPolynomial(coeffs=[0.963126,0.0601935,-2.65051e-05,4.89906e-09,-3.3308e-13,44456.9,36.0981], Tmin=(837.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Tertalkyl) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C([CH2])[CH]CCC(6210)',
    structure = SMILES('[CH2]C([CH2])[CH]CCC'),
    E0 = (383.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595958,0.072657,-4.89212e-05,1.88569e-08,-3.16737e-12,46232.6,33.699], Tmin=(100,'K'), Tmax=(1333.65,'K')), NASAPolynomial(coeffs=[9.97376,0.0445302,-1.72861e-05,3.04304e-09,-2.02971e-13,43731.3,-14.2435], Tmin=(1333.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])CC[CH]C(6237)',
    structure = SMILES('[CH2]C([CH2])CC[CH]C'),
    E0 = (383.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690633,0.0730084,-5.31014e-05,2.42699e-08,-5.00434e-12,46215.9,33.5436], Tmin=(100,'K'), Tmax=(1097.16,'K')), NASAPolynomial(coeffs=[7.62645,0.0477226,-1.85323e-05,3.26526e-09,-2.18313e-13,44693.9,-0.560973], Tmin=(1097.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)[CH][CH]CC(259)',
    structure = SMILES('[CH2]C(C)[CH][CH]CC'),
    E0 = (372.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,257.302,1053.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273575,'amu*angstrom^2'), symmetry=1, barrier=(0.964158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80979,0.0566356,1.2766e-05,-1.01655e-07,8.44898e-11,44901,30.9983], Tmin=(100,'K'), Tmax=(470.168,'K')), NASAPolynomial(coeffs=[3.19492,0.0566127,-2.46829e-05,4.64854e-09,-3.24377e-13,44640.7,23.9788], Tmin=(470.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2][C]([CH2])CCCC(6238)',
    structure = SMILES('[CH2][C]([CH2])CCCC'),
    E0 = (374.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521328,0.0808674,-8.27564e-05,6.07051e-08,-1.96934e-11,45132.8,32.5527], Tmin=(100,'K'), Tmax=(829.066,'K')), NASAPolynomial(coeffs=[5.35724,0.0522952,-2.15807e-05,3.88863e-09,-2.61753e-13,44511,11.2149], Tmin=(829.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCCC([CH2])[CH2](6239)',
    structure = SMILES('[CH2]CCCC([CH2])[CH2]'),
    E0 = (394.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,206.05,824.201,1632.34],'cm^-1')),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149345,'amu*angstrom^2'), symmetry=1, barrier=(3.57066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237229,0.0769853,-5.61277e-05,2.33219e-08,-4.0965e-12,47536.1,34.2199], Tmin=(100,'K'), Tmax=(1321.17,'K')), NASAPolynomial(coeffs=[12.2076,0.0407442,-1.49817e-05,2.55995e-09,-1.67867e-13,44373,-26.8646], Tmin=(1321.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C[CH][CH]C(260)',
    structure = SMILES('[CH2]C(C)C[CH][CH]C'),
    E0 = (372.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,225.632,1984.27,4000],'cm^-1')),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075197,'amu*angstrom^2'), symmetry=1, barrier=(2.35018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.963544,0.0713322,-6.29132e-05,4.61373e-08,-1.62613e-11,44924.2,34.0424], Tmin=(100,'K'), Tmax=(797.294,'K')), NASAPolynomial(coeffs=[2.76982,0.056446,-2.39495e-05,4.39515e-09,-2.99743e-13,44821.3,26.8982], Tmin=(797.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])C(262)',
    structure = SMILES('[CH2]C[CH]CC([CH2])C'),
    E0 = (383.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,217.504,512.993,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0041447,'amu*angstrom^2'), symmetry=1, barrier=(0.142352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995406,0.0695583,-4.56127e-05,1.82633e-08,-3.46478e-12,46223.3,32.9793], Tmin=(100,'K'), Tmax=(1119.2,'K')), NASAPolynomial(coeffs=[6.00483,0.0516547,-2.16175e-05,3.97022e-09,-2.72083e-13,45102,8.24756], Tmin=(1119.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(Isobutyl) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH]CC(6240)',
    structure = SMILES('[CH2][C]([CH2])C[CH]CC'),
    E0 = (568.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300.705,2672.4,3846.6],'cm^-1')),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735937,0.081846,-0.000108643,9.96047e-08,-3.62672e-11,68507.5,34.3431], Tmin=(100,'K'), Tmax=(876.764,'K')), NASAPolynomial(coeffs=[1.17354,0.056099,-2.39605e-05,4.31817e-09,-2.87443e-13,69343.6,37.4954], Tmin=(876.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJCC) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH]CC(6241)',
    structure = SMILES('[CH2]C([CH2])[CH][CH]CC'),
    E0 = (577.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350.462,3118.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023177,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10904,0.0522113,4.37709e-05,-2.09557e-07,1.92163e-10,69553.1,31.1515], Tmin=(100,'K'), Tmax=(421.018,'K')), NASAPolynomial(coeffs=[4.32379,0.0507898,-2.10676e-05,3.80124e-09,-2.55772e-13,69192.7,20.3174], Tmin=(421.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[CH2](6131)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH2]'),
    E0 = (612.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1171.04,1177.35],'cm^-1')),
        HinderedRotor(inertia=(0.00401418,'amu*angstrom^2'), symmetry=1, barrier=(3.85633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00405621,'amu*angstrom^2'), symmetry=1, barrier=(3.95424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171241,'amu*angstrom^2'), symmetry=1, barrier=(3.93717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16792,'amu*angstrom^2'), symmetry=1, barrier=(3.86082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169879,'amu*angstrom^2'), symmetry=1, barrier=(3.90586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64318,0.0551015,-3.22854e-05,-1.62113e-09,1.00768e-11,73724.5,29.706], Tmin=(100,'K'), Tmax=(647.81,'K')), NASAPolynomial(coeffs=[5.92381,0.0380635,-1.45841e-05,2.54567e-09,-1.69303e-13,72972.8,9.39165], Tmin=(647.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C[CH][CH]C(6242)',
    structure = SMILES('[CH2]C([CH2])C[CH][CH]C'),
    E0 = (577.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,293.508,1838.13,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00293027,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.87627,0.0742874,-7.99072e-05,6.43259e-08,-2.21288e-11,69592,35.4413], Tmin=(100,'K'), Tmax=(870.3,'K')), NASAPolynomial(coeffs=[3.13781,0.05207,-2.12368e-05,3.77352e-09,-2.50661e-13,69646.1,27.4173], Tmin=(870.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])[CH2](6243)',
    structure = SMILES('[CH2]C[CH]CC([CH2])[CH2]'),
    E0 = (588.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,236.79,965.874,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27905,0.0665007,-3.16217e-05,-2.44808e-08,3.07083e-11,70875.6,33.1397], Tmin=(100,'K'), Tmax=(566.504,'K')), NASAPolynomial(coeffs=[5.74833,0.0484593,-1.96383e-05,3.53141e-09,-2.38781e-13,70152.3,12.2036], Tmin=(566.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJCC) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C)C[CH]CC(247)',
    structure = SMILES('C=C(C)C[CH]CC'),
    E0 = (91.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.251,0.0650146,-3.38669e-05,8.35878e-09,-8.35867e-13,11076.5,27.7197], Tmin=(100,'K'), Tmax=(2093.64,'K')), NASAPolynomial(coeffs=[13.5653,0.0414877,-1.7011e-05,2.99144e-09,-1.94959e-13,5920.13,-40.7886], Tmin=(2093.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(=C)CCCC(6244)',
    structure = SMILES('[CH2]C(=C)CCCC'),
    E0 = (48.3619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291066,0.0699195,-2.75167e-05,-8.23717e-09,6.58711e-12,5959.97,28.2282], Tmin=(100,'K'), Tmax=(1084.29,'K')), NASAPolynomial(coeffs=[14.2315,0.0400031,-1.58886e-05,2.91008e-09,-2.01681e-13,1672.39,-45.9856], Tmin=(1084.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.3619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C=CCC(248)',
    structure = SMILES('[CH2]C(C)C=CCC'),
    E0 = (98.1373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,567.625],'cm^-1')),
        HinderedRotor(inertia=(0.121443,'amu*angstrom^2'), symmetry=1, barrier=(2.79221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121441,'amu*angstrom^2'), symmetry=1, barrier=(2.79217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121441,'amu*angstrom^2'), symmetry=1, barrier=(2.79218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616538,'amu*angstrom^2'), symmetry=1, barrier=(14.1754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619995,'amu*angstrom^2'), symmetry=1, barrier=(14.1754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619519,0.0628016,-1.18176e-05,-2.36944e-08,1.23304e-11,11935,30.7111], Tmin=(100,'K'), Tmax=(1012.26,'K')), NASAPolynomial(coeffs=[12.69,0.040009,-1.4947e-05,2.67144e-09,-1.83895e-13,8215.32,-33.9715], Tmin=(1012.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.1373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)CC=CC(249)',
    structure = SMILES('[CH2]C(C)CC=CC'),
    E0 = (95.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,247.363,247.366],'cm^-1')),
        HinderedRotor(inertia=(0.139652,'amu*angstrom^2'), symmetry=1, barrier=(6.06382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00275501,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139645,'amu*angstrom^2'), symmetry=1, barrier=(6.06382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00275511,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393575,'amu*angstrom^2'), symmetry=1, barrier=(17.0897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500022,0.0673883,-2.93978e-05,-2.14166e-09,3.90934e-12,11672.1,30.4409], Tmin=(100,'K'), Tmax=(1102.6,'K')), NASAPolynomial(coeffs=[12.0041,0.0415868,-1.59725e-05,2.84672e-09,-1.93297e-13,8166.76,-30.5748], Tmin=(1102.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C[CH]C(153)',
    structure = SMILES('[CH2]C([CH2])C[CH]C'),
    E0 = (407.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1161.67],'cm^-1')),
        HinderedRotor(inertia=(0.16627,'amu*angstrom^2'), symmetry=1, barrier=(3.82288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00399131,'amu*angstrom^2'), symmetry=1, barrier=(3.82163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166204,'amu*angstrom^2'), symmetry=1, barrier=(3.82136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00399205,'amu*angstrom^2'), symmetry=1, barrier=(3.82207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166113,'amu*angstrom^2'), symmetry=1, barrier=(3.81926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3039,0.0586266,-4.0687e-05,1.80742e-08,-3.65726e-12,49054.6,29.096], Tmin=(100,'K'), Tmax=(1112.4,'K')), NASAPolynomial(coeffs=[6.6473,0.039413,-1.47791e-05,2.54772e-09,-1.67926e-13,47865.8,2.74792], Tmin=(1112.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'CH2(S)(14)',
    structure = SMILES('[CH2]'),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144068,5.45069e-06,-3.58002e-09,7.56192e-13,50400.6,-0.411765], Tmin=(100,'K'), Tmax=(1442.36,'K')), NASAPolynomial(coeffs=[2.62648,0.00394763,-1.49924e-06,2.54539e-10,-1.62956e-14,50691.8,6.78378], Tmin=(1442.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][CH]CC[CH]CC(299)',
    structure = SMILES('[CH2][CH]CC[CH]CC'),
    E0 = (378.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,405.986,1418.08,2619.56,3831.92],'cm^-1')),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522728,'amu*angstrom^2'), symmetry=1, barrier=(5.35927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584014,-2.57862e-05,4.71004e-09,-3.04691e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.47,'K')), NASAPolynomial(coeffs=[41.6474,0.010469,-4.46141e-06,6.56552e-10,-3.30716e-14,20282.4,-206.655], Tmin=(2817.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])CC(239)',
    structure = SMILES('[CH2]C([CH2])C([CH2])CC'),
    E0 = (391.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1600],'cm^-1')),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3395.46,'J/mol'), sigma=(6.43099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.36 K, Pc=28.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354184,0.0739518,-4.30025e-05,4.70431e-09,4.28316e-12,47230.8,33.5164], Tmin=(100,'K'), Tmax=(923.294,'K')), NASAPolynomial(coeffs=[11.9522,0.0397531,-1.35139e-05,2.23689e-09,-1.45927e-13,44405.2,-25.216], Tmin=(923.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])CC([CH2])C(6245)',
    structure = SMILES('[CH2]C([CH2])CC([CH2])C'),
    E0 = (388.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1600],'cm^-1')),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161725,'amu*angstrom^2'), symmetry=1, barrier=(3.71837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354184,0.0739518,-4.30025e-05,4.70431e-09,4.28316e-12,46828.3,33.5164], Tmin=(100,'K'), Tmax=(923.294,'K')), NASAPolynomial(coeffs=[11.9522,0.0397531,-1.35139e-05,2.23689e-09,-1.45927e-13,44002.6,-25.216], Tmin=(923.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'CC[CH]CC1CC1(6246)',
    structure = SMILES('CC[CH]CC1CC1'),
    E0 = (132.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08625,0.0516162,1.18437e-05,-4.08156e-08,1.62853e-11,16035.8,28.6372], Tmin=(100,'K'), Tmax=(1051.77,'K')), NASAPolynomial(coeffs=[10.4395,0.0441488,-1.75878e-05,3.24532e-09,-2.26577e-13,12513.8,-24.3485], Tmin=(1051.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C1CC(CC)C1(6086)',
    structure = SMILES('[CH2]C1CC(CC)C1'),
    E0 = (129.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21259,0.0406227,6.16011e-05,-1.04592e-07,4.19798e-11,15698,26.4384], Tmin=(100,'K'), Tmax=(962.846,'K')), NASAPolynomial(coeffs=[14.4625,0.036959,-1.27368e-05,2.30192e-09,-1.65181e-13,10764.8,-49.3509], Tmin=(962.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C[CH]CC(220)',
    structure = SMILES('[CH2][CH]C[CH]CC'),
    E0 = (402.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,214.497,397.673,3998.61],'cm^-1')),
        HinderedRotor(inertia=(0.00686178,'amu*angstrom^2'), symmetry=1, barrier=(15.0313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00686178,'amu*angstrom^2'), symmetry=1, barrier=(15.0313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00686178,'amu*angstrom^2'), symmetry=1, barrier=(15.0313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00686178,'amu*angstrom^2'), symmetry=1, barrier=(15.0313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00686178,'amu*angstrom^2'), symmetry=1, barrier=(15.0313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62267,0.0581169,-5.74867e-05,5.14491e-08,-2.07613e-11,48463.8,29.7123], Tmin=(100,'K'), Tmax=(801.747,'K')), NASAPolynomial(coeffs=[0.541337,0.0516173,-2.30726e-05,4.32877e-09,-2.98228e-13,49019.5,37.0742], Tmin=(801.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH]CC(45)',
    structure = SMILES('[CH]CC'),
    E0 = (327.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,319.066,319.142,1779.44],'cm^-1')),
        HinderedRotor(inertia=(0.00165544,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165526,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00517,0.0155111,1.85851e-05,-3.16862e-08,1.23732e-11,39453.6,11.9344], Tmin=(100,'K'), Tmax=(982.292,'K')), NASAPolynomial(coeffs=[6.73204,0.0159276,-5.86166e-06,1.06538e-09,-7.51285e-14,37969.2,-9.8082], Tmin=(982.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])[CH2](6116)',
    structure = SMILES('[CH2]C([CH2])[CH2]'),
    E0 = (461.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00529301,'amu*angstrom^2'), symmetry=1, barrier=(3.08377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00523118,'amu*angstrom^2'), symmetry=1, barrier=(3.03916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130467,'amu*angstrom^2'), symmetry=1, barrier=(2.9997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9185,0.0363463,-2.61059e-05,1.15462e-08,-2.04218e-12,55637,20.3501], Tmin=(100,'K'), Tmax=(1631.3,'K')), NASAPolynomial(coeffs=[7.087,0.0180707,-4.14987e-06,4.68151e-10,-2.1815e-14,54696.1,-4.82954], Tmin=(1631.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH2](6)',
    structure = SMILES('C[CH2]'),
    E0 = (108.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,474.132,1048.55,2319.88,2320.55,2321.73],'cm^-1')),
        HinderedRotor(inertia=(0.000749852,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82183,-0.00343357,5.09256e-05,-6.2021e-08,2.37073e-11,13066,7.61644], Tmin=(100,'K'), Tmax=(900.314,'K')), NASAPolynomial(coeffs=[5.15622,0.00943121,-1.81945e-06,2.21194e-10,-1.4348e-14,12064.1,-2.91103], Tmin=(900.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ)"""),
)

species(
    label = '[CH]CC([CH2])[CH2](6143)',
    structure = SMILES('[CH]CC([CH2])[CH2]'),
    E0 = (684.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,994.78,996.827,996.949,999.352],'cm^-1')),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.52207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106583,'amu*angstrom^2'), symmetry=1, barrier=(2.45054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104897,'amu*angstrom^2'), symmetry=1, barrier=(2.4118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107425,'amu*angstrom^2'), symmetry=1, barrier=(2.46992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35611,0.0506983,-4.07186e-05,1.87001e-08,-3.49068e-12,82440,24.9335], Tmin=(100,'K'), Tmax=(1370.02,'K')), NASAPolynomial(coeffs=[10.6767,0.0214606,-6.49018e-06,9.65436e-10,-5.76292e-14,80076.1,-22.2738], Tmin=(1370.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C[C]CC(6247)',
    structure = SMILES('[CH2]C([CH2])C[C]CC'),
    E0 = (637.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.172908,0.0790275,-6.43323e-05,2.98526e-08,-5.767e-12,76760.4,32.5689], Tmin=(100,'K'), Tmax=(1224.8,'K')), NASAPolynomial(coeffs=[12.8027,0.0377806,-1.38176e-05,2.35704e-09,-1.5472e-13,73666.6,-30.9233], Tmin=(1224.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C([CH2])C[CH]CC(6248)',
    structure = SMILES('[CH]C([CH2])C[CH]CC'),
    E0 = (626.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,215.054,814.885,1018.61,1222.33,1504.19,1719.08],'cm^-1')),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138907,'amu*angstrom^2'), symmetry=1, barrier=(3.52051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839062,0.0725934,-5.68797e-05,2.82419e-08,-6.40776e-12,75451.7,32.7633], Tmin=(100,'K'), Tmax=(996.49,'K')), NASAPolynomial(coeffs=[6.9056,0.0482417,-2.02232e-05,3.71798e-09,-2.55158e-13,74242.7,3.51716], Tmin=(996.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    E0 = (383.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (465.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (520.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (518.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (511.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (433.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (482.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (495.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (524.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (542.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (542.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (501.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (527.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (535.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (429.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (436.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (739.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (780.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (789.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (747.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (789.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (800.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (406.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (446.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (446.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (408.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (826.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (540.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (551.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (582.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (391.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (391.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (817.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (823.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (827.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (848.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (838.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['C=CCC(36)', '[CH2]C=C(87)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2]C(=C)C[CH]CC(6234)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C([CH2])C=CCC(6235)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C([CH2])CC=CC(6236)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(28)', 'C=CC[CH]CC(216)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C(87)', '[CH2][CH]CC(39)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00179511,'m^3/(mol*s)'), n=2.50446, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CCC(36)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0172287,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH3](11)', '[CH2]C([CH2])CC=C(6127)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2][C](C)C[CH]CC(258)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C([CH2])[CH]CCC(6210)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C([CH2])CC[CH]C(6237)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C(C)[CH][CH]CC(259)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2][C]([CH2])CCCC(6238)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CCCC([CH2])[CH2](6239)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C(C)C[CH][CH]C(260)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(524000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C[CH]CC([CH2])C(262)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(10873.3,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R6HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC(39)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2][C]([CH2])C[CH]CC(6240)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2]C([CH2])[CH][CH]CC(6241)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH3](11)', '[CH2][CH]CC([CH2])[CH2](6131)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C([CH2])C[CH][CH]C(6242)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2]C[CH]CC([CH2])[CH2](6243)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['C=C(C)C[CH]CC(247)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C(=C)CCCC(6244)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C(C)C=CCC(248)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C(C)CC=CC(249)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.50442e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])C[CH]C(153)', 'CH2(S)(14)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2][CH]CC[CH]CC(299)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])C([CH2])CC(239)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH2])CC([CH2])C(6245)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['CC[CH]CC1CC1(6246)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH2])C[CH]CC(261)'],
    products = ['[CH2]C1CC(CC)C1(6086)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C[CH]CC(220)', 'CH2(T)(28)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]CC(45)', '[CH2]C([CH2])[CH2](6116)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.44562e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH2](6)', '[CH]CC([CH2])[CH2](6143)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', '[CH2]C([CH2])C[C]CC(6247)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', '[CH]C([CH2])C[CH]CC(6248)'],
    products = ['[CH2]C([CH2])C[CH]CC(261)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '1728',
    isomers = [
        '[CH2]C([CH2])C[CH]CC(261)',
    ],
    reactants = [
        ('C=CCC(36)', '[CH2]C=C(87)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '1728',
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


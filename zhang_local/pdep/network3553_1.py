species(
    label = '[CH2]C(O)C([CH2])C=O(12792)',
    structure = SMILES('[CH2]C(O)C([CH2])C=O'),
    E0 = (-7.01478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,238.306],'cm^-1')),
        HinderedRotor(inertia=(0.00296845,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00296846,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206246,'amu*angstrom^2'), symmetry=1, barrier=(8.31156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206246,'amu*angstrom^2'), symmetry=1, barrier=(8.31156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347604,'amu*angstrom^2'), symmetry=1, barrier=(14.0082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428534,0.0833066,-0.000109453,8.16414e-08,-2.48125e-11,-719.471,29.9127], Tmin=(100,'K'), Tmax=(800.968,'K')), NASAPolynomial(coeffs=[10.5636,0.0326932,-1.46686e-05,2.75145e-09,-1.89506e-13,-2343.07,-16.734], Tmin=(800.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.01478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CJCO)"""),
)

species(
    label = 'C=CO(576)',
    structure = SMILES('C=CO'),
    E0 = (-166.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.24798,'amu*angstrom^2'), symmetry=1, barrier=(28.6936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92544,0.00836062,5.0344e-05,-8.45232e-08,3.72335e-11,-19989.5,9.0776], Tmin=(100,'K'), Tmax=(898.452,'K')), NASAPolynomial(coeffs=[15.1116,-0.00538397,5.65903e-06,-1.18193e-09,7.91212e-14,-23814.2,-57.5076], Tmin=(898.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)C1CC1[O](13855)',
    structure = SMILES('[CH2]C(O)C1CC1[O]'),
    E0 = (82.4797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835562,0.05613,-1.34031e-05,-2.8538e-08,1.65533e-11,10046.1,26.9445], Tmin=(100,'K'), Tmax=(964.369,'K')), NASAPolynomial(coeffs=[17.3612,0.0196435,-6.51567e-06,1.172e-09,-8.48853e-14,5367.99,-59.9117], Tmin=(964.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.4797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1C([O])CC1O(13856)',
    structure = SMILES('[CH2]C1C([O])CC1O'),
    E0 = (71.5536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30698,0.0434367,2.07792e-05,-6.19843e-08,2.81882e-11,8717.26,25.7218], Tmin=(100,'K'), Tmax=(946.659,'K')), NASAPolynomial(coeffs=[15.8862,0.0212168,-6.41589e-06,1.11342e-09,-8.06778e-14,4192.28,-53.1357], Tmin=(946.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.5536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C(O)C(=C)C=O(13857)',
    structure = SMILES('[CH2]C(O)C(=C)C=O'),
    E0 = (-108.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,348.711],'cm^-1')),
        HinderedRotor(inertia=(0.155005,'amu*angstrom^2'), symmetry=1, barrier=(13.3752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00138607,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154994,'amu*angstrom^2'), symmetry=1, barrier=(13.3741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154895,'amu*angstrom^2'), symmetry=1, barrier=(13.3747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16011,0.0633355,-5.4549e-05,2.38375e-08,-4.27671e-12,-12953.7,26.7966], Tmin=(100,'K'), Tmax=(1298,'K')), NASAPolynomial(coeffs=[12.7767,0.0275371,-1.31796e-05,2.58981e-09,-1.84321e-13,-15969.4,-32.2767], Tmin=(1298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(C=O)C(=C)O(13858)',
    structure = SMILES('[CH2]C(C=O)C(=C)O'),
    E0 = (-134.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,272.621],'cm^-1')),
        HinderedRotor(inertia=(0.244216,'amu*angstrom^2'), symmetry=1, barrier=(12.8679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244204,'amu*angstrom^2'), symmetry=1, barrier=(12.8674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24407,'amu*angstrom^2'), symmetry=1, barrier=(12.866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244543,'amu*angstrom^2'), symmetry=1, barrier=(12.8685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.513866,0.0757074,-8.45025e-05,4.91838e-08,-1.13575e-11,-15995.4,26.4734], Tmin=(100,'K'), Tmax=(1056.23,'K')), NASAPolynomial(coeffs=[14.5032,0.0227274,-9.26068e-06,1.69152e-09,-1.16207e-13,-18950.5,-41.782], Tmin=(1056.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][CH]O(578)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (135.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0943891,'amu*angstrom^2'), symmetry=1, barrier=(7.36374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0943883,'amu*angstrom^2'), symmetry=1, barrier=(7.36374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67796,0.0299043,-3.92791e-05,2.86662e-08,-8.27893e-12,16321.6,12.7466], Tmin=(100,'K'), Tmax=(918.072,'K')), NASAPolynomial(coeffs=[6.82026,0.00999271,-3.7014e-06,6.19844e-10,-3.95112e-14,15639.5,-6.45576], Tmin=(918.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CJCO)"""),
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
    label = '[CH2]C(O)C=C(737)',
    structure = SMILES('[CH2]C(O)C=C'),
    E0 = (22.2606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,315.491],'cm^-1')),
        HinderedRotor(inertia=(0.194591,'amu*angstrom^2'), symmetry=1, barrier=(13.7996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00167931,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192492,'amu*angstrom^2'), symmetry=1, barrier=(13.7896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79026,0.0401682,-1.22391e-05,-1.4209e-08,8.68931e-12,2764.37,21.9101], Tmin=(100,'K'), Tmax=(989.738,'K')), NASAPolynomial(coeffs=[12.1514,0.017266,-6.28282e-06,1.14653e-09,-8.14646e-14,-215.826,-32.6638], Tmin=(989.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO)"""),
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
    label = '[CH2]C(C=C)C=O(12631)',
    structure = SMILES('[CH2]C(C=C)C=O'),
    E0 = (80.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,276.463],'cm^-1')),
        HinderedRotor(inertia=(0.138121,'amu*angstrom^2'), symmetry=1, barrier=(7.50399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140152,'amu*angstrom^2'), symmetry=1, barrier=(7.5006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233274,'amu*angstrom^2'), symmetry=1, barrier=(12.7234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6842,0.054824,-5.37097e-05,3.28464e-08,-8.87785e-12,9763.62,22.2265], Tmin=(100,'K'), Tmax=(864.442,'K')), NASAPolynomial(coeffs=[6.53096,0.0323977,-1.47964e-05,2.83722e-09,-1.99414e-13,8925.64,-0.45052], Tmin=(864.442,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(O)C(C)=C[O](13859)',
    structure = SMILES('[CH2]C(O)C(C)=C[O]'),
    E0 = (-84.1244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,218.189,221.693],'cm^-1')),
        HinderedRotor(inertia=(0.00345228,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501507,'amu*angstrom^2'), symmetry=1, barrier=(17.1371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.490554,'amu*angstrom^2'), symmetry=1, barrier=(17.1658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488734,'amu*angstrom^2'), symmetry=1, barrier=(17.1218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.33405,0.0688206,-4.34862e-05,-2.44425e-09,8.58458e-12,-9975.17,29.1624], Tmin=(100,'K'), Tmax=(960.401,'K')), NASAPolynomial(coeffs=[19.4233,0.0168553,-5.33778e-06,9.32846e-10,-6.67768e-14,-14911.9,-68.7733], Tmin=(960.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.1244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(C=O)[C](C)O(13684)',
    structure = SMILES('[CH2]C(C=O)[C](C)O'),
    E0 = (-41.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332726,0.0880793,-0.000131657,1.12039e-07,-3.78706e-11,-4923.62,29.0395], Tmin=(100,'K'), Tmax=(840.337,'K')), NASAPolynomial(coeffs=[8.64145,0.0360691,-1.6577e-05,3.09678e-09,-2.10707e-13,-5880.08,-6.98215], Tmin=(840.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][C](O)C(C)C=O(13860)',
    structure = SMILES('[CH2][C](O)C(C)C=O'),
    E0 = (-40.8977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.470784,0.0831867,-0.000115373,9.30552e-08,-3.04797e-11,-4797.08,28.8281], Tmin=(100,'K'), Tmax=(813.811,'K')), NASAPolynomial(coeffs=[9.20671,0.0346627,-1.56387e-05,2.92011e-09,-1.9965e-13,-6033.99,-10.3814], Tmin=(813.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.8977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C(C)[C]=O(13861)',
    structure = SMILES('[CH2]C(O)C(C)[C]=O'),
    E0 = (-58.8369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.832585,0.0712633,-7.32232e-05,4.17783e-08,-9.8121e-12,-6963.69,29.2671], Tmin=(100,'K'), Tmax=(1020.02,'K')), NASAPolynomial(coeffs=[11.2722,0.0303249,-1.30217e-05,2.43225e-09,-1.68764e-13,-9093.44,-21.305], Tmin=(1020.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.8369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(C=O)C(C)[O](12750)',
    structure = SMILES('[CH2]C(C=O)C(C)[O]'),
    E0 = (11.7569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,250.745,250.875],'cm^-1')),
        HinderedRotor(inertia=(0.00269962,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00269174,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25729,'amu*angstrom^2'), symmetry=1, barrier=(11.4093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254878,'amu*angstrom^2'), symmetry=1, barrier=(11.4127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4030.69,'J/mol'), sigma=(6.74566,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.58 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907792,0.0720923,-7.6221e-05,4.68729e-08,-1.21665e-11,1521.9,27.3735], Tmin=(100,'K'), Tmax=(916.66,'K')), NASAPolynomial(coeffs=[9.34716,0.0352656,-1.59586e-05,3.04529e-09,-2.13453e-13,-25.3057,-12.6072], Tmin=(916.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.7569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(=C[O])C(C)O(13687)',
    structure = SMILES('[CH2]C(=C[O])C(C)O'),
    E0 = (-144.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499589,0.0609652,-1.45013e-05,-3.60057e-08,2.12808e-11,-17204.3,28.2155], Tmin=(100,'K'), Tmax=(946.941,'K')), NASAPolynomial(coeffs=[20.8059,0.0147468,-3.95132e-06,6.82384e-10,-5.198e-14,-22823.7,-78.0087], Tmin=(946.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)C(C)O(13688)',
    structure = SMILES('[CH2]C([C]=O)C(C)O'),
    E0 = (-59.9152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728736,0.0756568,-8.73348e-05,5.73964e-08,-1.55558e-11,-7091.57,29.3623], Tmin=(100,'K'), Tmax=(889.8,'K')), NASAPolynomial(coeffs=[10.3606,0.0323557,-1.43355e-05,2.70037e-09,-1.87583e-13,-8805.58,-15.9809], Tmin=(889.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.9152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([O])C(C)C=O(13689)',
    structure = SMILES('[CH2]C([O])C(C)C=O'),
    E0 = (12.8353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0133696,'amu*angstrom^2'), symmetry=1, barrier=(2.10586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133493,'amu*angstrom^2'), symmetry=1, barrier=(2.06266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761204,'amu*angstrom^2'), symmetry=1, barrier=(11.6816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507167,'amu*angstrom^2'), symmetry=1, barrier=(11.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996761,0.0679205,-6.30765e-05,3.27318e-08,-7.13073e-12,1650.32,27.3283], Tmin=(100,'K'), Tmax=(1080.86,'K')), NASAPolynomial(coeffs=[10.4702,0.0328623,-1.44243e-05,2.72406e-09,-1.90165e-13,-397.611,-19.1124], Tmin=(1080.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.8353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C=O(12931)',
    structure = SMILES('[CH2][CH]C([CH2])C=O'),
    E0 = (359.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,567.782],'cm^-1')),
        HinderedRotor(inertia=(0.000287703,'amu*angstrom^2'), symmetry=1, barrier=(2.0696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895652,'amu*angstrom^2'), symmetry=1, barrier=(2.05928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173514,'amu*angstrom^2'), symmetry=1, barrier=(3.98942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54319,'amu*angstrom^2'), symmetry=1, barrier=(35.4809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65561,0.0546676,-5.19115e-05,2.93341e-08,-7.17287e-12,43261.3,26.8024], Tmin=(100,'K'), Tmax=(956.996,'K')), NASAPolynomial(coeffs=[7.52151,0.0301497,-1.34821e-05,2.56333e-09,-1.79451e-13,42138.6,-1.23933], Tmin=(956.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CJC(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C=O(12644)',
    structure = SMILES('[CH2]C([O])C([CH2])C=O'),
    E0 = (223.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,237.377,2887.88],'cm^-1')),
        HinderedRotor(inertia=(0.31931,'amu*angstrom^2'), symmetry=1, barrier=(12.7646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00299155,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319235,'amu*angstrom^2'), symmetry=1, barrier=(12.7648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788709,'amu*angstrom^2'), symmetry=1, barrier=(31.5423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624896,0.0788821,-0.000101868,7.41622e-08,-2.20348e-11,26979.6,29.0184], Tmin=(100,'K'), Tmax=(817.816,'K')), NASAPolynomial(coeffs=[10.399,0.0310773,-1.41882e-05,2.68947e-09,-1.86674e-13,25380.9,-16.1705], Tmin=(817.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)OJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])O(729)',
    structure = SMILES('[CH2][CH]C([CH2])O'),
    E0 = (299.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,944.581],'cm^-1')),
        HinderedRotor(inertia=(0.17604,'amu*angstrom^2'), symmetry=1, barrier=(4.04751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00639217,'amu*angstrom^2'), symmetry=1, barrier=(4.04713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00639262,'amu*angstrom^2'), symmetry=1, barrier=(4.04713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.544563,'amu*angstrom^2'), symmetry=1, barrier=(12.5206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55872,0.0478807,-3.98632e-05,1.75683e-08,-3.11477e-12,36170,25.6447], Tmin=(100,'K'), Tmax=(1352.88,'K')), NASAPolynomial(coeffs=[11.7473,0.0177563,-6.46257e-06,1.10914e-09,-7.3245e-14,33413.2,-26.5888], Tmin=(1352.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(=C[O])C([CH2])O(13862)',
    structure = SMILES('[CH2]C(=C[O])C([CH2])O'),
    E0 = (67.3747,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,457.843,457.895],'cm^-1')),
        HinderedRotor(inertia=(0.129876,'amu*angstrom^2'), symmetry=1, barrier=(19.3228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129725,'amu*angstrom^2'), symmetry=1, barrier=(19.3301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129674,'amu*angstrom^2'), symmetry=1, barrier=(19.3277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130071,'amu*angstrom^2'), symmetry=1, barrier=(19.3261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.350386,0.0661182,-3.41272e-05,-1.70005e-08,1.51869e-11,8247.68,29.3848], Tmin=(100,'K'), Tmax=(940.361,'K')), NASAPolynomial(coeffs=[21.5872,0.011021,-2.44831e-06,3.89784e-10,-3.04461e-14,2695.63,-80.0488], Tmin=(940.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.3747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CJCO) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C](O)C([CH2])C=O(13863)',
    structure = SMILES('[CH2][C](O)C([CH2])C=O'),
    E0 = (169.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,360,370,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,257.359],'cm^-1')),
        HinderedRotor(inertia=(0.00254521,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152489,'amu*angstrom^2'), symmetry=1, barrier=(7.16711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15249,'amu*angstrom^2'), symmetry=1, barrier=(7.16711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152488,'amu*angstrom^2'), symmetry=1, barrier=(7.16711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152487,'amu*angstrom^2'), symmetry=1, barrier=(7.16711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15077,0.0936222,-0.000152638,1.32727e-07,-4.46161e-11,20529.8,30.3259], Tmin=(100,'K'), Tmax=(865.691,'K')), NASAPolynomial(coeffs=[9.60845,0.032022,-1.48864e-05,2.75946e-09,-1.85438e-13,19563,-10.0639], Tmin=(865.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C([CH2])[C]=O(13864)',
    structure = SMILES('[CH2]C(O)C([CH2])[C]=O'),
    E0 = (151.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,4000],'cm^-1')),
        HinderedRotor(inertia=(0.11571,'amu*angstrom^2'), symmetry=1, barrier=(9.29407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204647,'amu*angstrom^2'), symmetry=1, barrier=(16.4503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433915,'amu*angstrom^2'), symmetry=1, barrier=(34.8385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115482,'amu*angstrom^2'), symmetry=1, barrier=(9.29484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11552,'amu*angstrom^2'), symmetry=1, barrier=(9.29518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448961,0.0824386,-0.000113091,8.50205e-08,-2.56498e-11,18366,30.9942], Tmin=(100,'K'), Tmax=(811.147,'K')), NASAPolynomial(coeffs=[11.4779,0.028044,-1.24886e-05,2.32552e-09,-1.59163e-13,16577,-19.9041], Tmin=(811.147,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CC(C)CJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(O)C1[CH]OC1(13865)',
    structure = SMILES('[CH2]C(O)C1[CH]OC1'),
    E0 = (66.6533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.332447,0.0739207,-6.95383e-05,3.42484e-08,-6.31357e-12,8190.9,29.6712], Tmin=(100,'K'), Tmax=(1593.05,'K')), NASAPolynomial(coeffs=[16.0138,0.017551,-2.03048e-06,-4.15844e-11,1.5326e-14,4927.52,-50.6981], Tmin=(1593.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.6533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C1[CH]OCC1O(13866)',
    structure = SMILES('[CH2]C1[CH]OCC1O'),
    E0 = (-18.3159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31368,0.0451882,2.08746e-05,-7.36416e-08,3.80555e-11,-2092.72,22.0257], Tmin=(100,'K'), Tmax=(855.487,'K')), NASAPolynomial(coeffs=[16.5721,0.0154801,-3.91043e-08,-4.53614e-10,4.25431e-14,-6226.96,-58.1104], Tmin=(855.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.3159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C=O)C(C)O(13694)',
    structure = SMILES('C=C(C=O)C(C)O'),
    E0 = (-320.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06664,0.0609153,-4.37886e-05,1.52398e-08,-2.13785e-12,-38395.2,26.505], Tmin=(100,'K'), Tmax=(1638.11,'K')), NASAPolynomial(coeffs=[15.3432,0.0260541,-1.18664e-05,2.24819e-09,-1.55123e-13,-43072.5,-49.4171], Tmin=(1638.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C(O)C(C)C=O(13867)',
    structure = SMILES('C=C(O)C(C)C=O'),
    E0 = (-344.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460395,0.0696695,-6.2508e-05,2.90112e-08,-5.34712e-12,-41306.1,26.315], Tmin=(100,'K'), Tmax=(1313.45,'K')), NASAPolynomial(coeffs=[16.0937,0.0220601,-8.13705e-06,1.41439e-09,-9.44371e-14,-45412.9,-53.3692], Tmin=(1313.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](C[O])C([CH2])O(9867)',
    structure = SMILES('[CH2][C](C[O])C([CH2])O'),
    E0 = (288.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,1380,1390,370,380,2900,435,257.302,908.501,2005.25],'cm^-1')),
        HinderedRotor(inertia=(0.0880245,'amu*angstrom^2'), symmetry=1, barrier=(3.34971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880245,'amu*angstrom^2'), symmetry=1, barrier=(3.34971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880245,'amu*angstrom^2'), symmetry=1, barrier=(3.34971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880245,'amu*angstrom^2'), symmetry=1, barrier=(3.34971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880245,'amu*angstrom^2'), symmetry=1, barrier=(3.34971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824883,0.0742635,-9.24954e-05,7.16916e-08,-2.31729e-11,34853.6,31.501], Tmin=(100,'K'), Tmax=(814.857,'K')), NASAPolynomial(coeffs=[7.64372,0.0367288,-1.59232e-05,2.92712e-09,-1.98858e-13,33877.2,0.827685], Tmin=(814.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C([CH2])C[O](13868)',
    structure = SMILES('[CH2][C](O)C([CH2])C[O]'),
    E0 = (312.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,1380,1390,370,380,2900,435,208.256,808.256,1666.05],'cm^-1')),
        HinderedRotor(inertia=(0.14745,'amu*angstrom^2'), symmetry=1, barrier=(3.58388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14745,'amu*angstrom^2'), symmetry=1, barrier=(3.58388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14745,'amu*angstrom^2'), symmetry=1, barrier=(3.58388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14745,'amu*angstrom^2'), symmetry=1, barrier=(3.58388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14745,'amu*angstrom^2'), symmetry=1, barrier=(3.58388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338779,0.0903953,-0.000143672,1.27435e-07,-4.35662e-11,37759,32.4722], Tmin=(100,'K'), Tmax=(881.008,'K')), NASAPolynomial(coeffs=[7.00962,0.038028,-1.69182e-05,3.07164e-09,-2.03833e-13,37440.5,5.99753], Tmin=(881.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C([CH2])[CH]O(13869)',
    structure = SMILES('[CH2][C](O)C([CH2])[CH]O'),
    E0 = (267.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.273147,0.101652,-0.00016177,1.33536e-07,-4.21398e-11,32322,33.3569], Tmin=(100,'K'), Tmax=(913.644,'K')), NASAPolynomial(coeffs=[12.223,0.0285781,-1.16496e-05,1.99636e-09,-1.26676e-13,30805,-21.6065], Tmin=(913.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCO) + radical(C2CsJOH) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C[O](13698)',
    structure = SMILES('[CH2]C([O])C([CH2])C[O]'),
    E0 = (366.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,990.419,990.451,990.48],'cm^-1')),
        HinderedRotor(inertia=(0.155487,'amu*angstrom^2'), symmetry=1, barrier=(3.57496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155361,'amu*angstrom^2'), symmetry=1, barrier=(3.57206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15546,'amu*angstrom^2'), symmetry=1, barrier=(3.57434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155513,'amu*angstrom^2'), symmetry=1, barrier=(3.57555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712346,0.0769615,-9.80807e-05,7.65611e-08,-2.47525e-11,44213.1,31.5178], Tmin=(100,'K'), Tmax=(816.54,'K')), NASAPolynomial(coeffs=[8.05093,0.0366322,-1.59493e-05,2.93585e-09,-1.99492e-13,43160.6,-1.50515], Tmin=(816.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])C([CH2])[CH]O(13699)',
    structure = SMILES('[CH2]C([O])C([CH2])[CH]O'),
    E0 = (321.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22839,0.0864709,-0.000108747,7.06188e-08,-1.68173e-11,38770.6,31.9575], Tmin=(100,'K'), Tmax=(743.957,'K')), NASAPolynomial(coeffs=[13.1897,0.0273313,-1.07761e-05,1.88477e-09,-1.2445e-13,36550.2,-28.7015], Tmin=(743.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCsJOH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[C-]#[O+](374)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (299.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33667,0.00896487,-2.66756e-05,3.61071e-08,-1.57199e-11,36069.2,-1.20266], Tmin=(100,'K'), Tmax=(865.594,'K')), NASAPolynomial(coeffs=[-0.394107,0.0117562,-6.47408e-06,1.26375e-09,-8.67562e-14,37256.3,19.3844], Tmin=(865.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.89,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH2]CC([CH2])O(3672)',
    structure = SMILES('[CH2]CC([CH2])O'),
    E0 = (100.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1380,1390,370,380,2900,435,702.232],'cm^-1')),
        HinderedRotor(inertia=(0.0749706,'amu*angstrom^2'), symmetry=1, barrier=(1.72372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0754573,'amu*angstrom^2'), symmetry=1, barrier=(1.73491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075095,'amu*angstrom^2'), symmetry=1, barrier=(1.72658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516648,'amu*angstrom^2'), symmetry=1, barrier=(11.8788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3512.21,'J/mol'), sigma=(6.17291,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.60 K, Pc=33.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53762,0.0510969,-4.16986e-05,1.84751e-08,-3.37147e-12,12125.7,23.2207], Tmin=(100,'K'), Tmax=(1292.36,'K')), NASAPolynomial(coeffs=[10.6393,0.0229262,-9.0017e-06,1.60832e-09,-1.08677e-13,9773.16,-23.024], Tmin=(1292.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)CC=C[O](12790)',
    structure = SMILES('[CH2]C(O)CC=C[O]'),
    E0 = (-68.5134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406561,0.06772,-4.23849e-05,-2.15442e-09,8.19593e-12,-8100.65,29.4809], Tmin=(100,'K'), Tmax=(960.095,'K')), NASAPolynomial(coeffs=[18.6882,0.0179115,-5.7463e-06,9.98562e-10,-7.06811e-14,-12825.9,-64.2993], Tmin=(960.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.5134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)[CH]CC=O(13870)',
    structure = SMILES('[CH2]C(O)[CH]CC=O'),
    E0 = (-11.8629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0362,0.0684816,-6.91899e-05,4.03026e-08,-9.90407e-12,-1322.82,30.0185], Tmin=(100,'K'), Tmax=(964.992,'K')), NASAPolynomial(coeffs=[9.45225,0.0335961,-1.49633e-05,2.83998e-09,-1.98645e-13,-2947.1,-10.2842], Tmin=(964.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.8629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(C=O)C[CH]O(12788)',
    structure = SMILES('[CH2]C(C=O)C[CH]O'),
    E0 = (-25.6631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.254097,0.0898527,-0.000135261,1.14743e-07,-3.85549e-11,-2958.83,29.792], Tmin=(100,'K'), Tmax=(843.473,'K')), NASAPolynomial(coeffs=[9.05553,0.035536,-1.62985e-05,3.03843e-09,-2.06355e-13,-3996.17,-8.5194], Tmin=(843.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.6631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCsJOH)"""),
)

species(
    label = 'O=CC1CCC1O(12794)',
    structure = SMILES('O=CC1CCC1O'),
    E0 = (-270.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53911,0.0408359,1.59778e-05,-4.68821e-08,1.98405e-11,-32486,23.8091], Tmin=(100,'K'), Tmax=(1001.8,'K')), NASAPolynomial(coeffs=[12.7864,0.0272467,-1.05693e-05,1.99083e-09,-1.43404e-13,-36311.1,-38.3169], Tmin=(1001.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C=COC([CH2])O(12791)',
    structure = SMILES('[CH2]C=COC([CH2])O'),
    E0 = (-58.6753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0808702,0.0707801,-3.58663e-05,-1.90901e-08,1.6557e-11,-6901.79,28.6299], Tmin=(100,'K'), Tmax=(942.434,'K')), NASAPolynomial(coeffs=[23.0054,0.0116932,-2.64171e-06,4.30212e-10,-3.39443e-14,-12919.7,-89.6115], Tmin=(942.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.6753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(=CO)C([CH2])O(13871)',
    structure = SMILES('[CH2]C(=CO)C([CH2])O'),
    E0 = (-74.0879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,350,440,435,1725,1380,1390,370,380,2900,435,294.489],'cm^-1')),
        HinderedRotor(inertia=(0.308833,'amu*angstrom^2'), symmetry=1, barrier=(18.9496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310782,'amu*angstrom^2'), symmetry=1, barrier=(18.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309175,'amu*angstrom^2'), symmetry=1, barrier=(18.9447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308305,'amu*angstrom^2'), symmetry=1, barrier=(18.9557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3088,'amu*angstrom^2'), symmetry=1, barrier=(18.9431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0594973,0.0730571,-3.73506e-05,-2.35822e-08,2.01603e-11,-8749.46,29.4467], Tmin=(100,'K'), Tmax=(918.956,'K')), NASAPolynomial(coeffs=[24.8859,0.00725829,2.18046e-07,-1.74995e-10,9.97478e-15,-15140.7,-98.6212], Tmin=(918.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.0879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CJCO) + radical(Allyl_P)"""),
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
    label = '[CH2]C(O)C=C[O](1110)',
    structure = SMILES('[CH2]C(O)C=C[O]'),
    E0 = (-45.0694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,344.234,345.179],'cm^-1')),
        HinderedRotor(inertia=(0.204319,'amu*angstrom^2'), symmetry=1, barrier=(17.2187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205606,'amu*angstrom^2'), symmetry=1, barrier=(17.2166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206412,'amu*angstrom^2'), symmetry=1, barrier=(17.1956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18004,0.0487819,-1.27556e-05,-2.98389e-08,1.8178e-11,-5306.9,25.3697], Tmin=(100,'K'), Tmax=(936.532,'K')), NASAPolynomial(coeffs=[18.491,0.00827652,-1.42509e-06,2.11294e-10,-1.83928e-14,-10015.5,-64.838], Tmin=(936.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.0694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]O)C=O(13541)',
    structure = SMILES('[CH2]C([CH]O)C=O'),
    E0 = (-1.88283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,324.475],'cm^-1')),
        HinderedRotor(inertia=(0.100253,'amu*angstrom^2'), symmetry=1, barrier=(7.47933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974971,'amu*angstrom^2'), symmetry=1, barrier=(7.48513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0991176,'amu*angstrom^2'), symmetry=1, barrier=(7.46516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980018,'amu*angstrom^2'), symmetry=1, barrier=(7.49914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932312,0.0746793,-0.000119922,1.04441e-07,-3.52704e-11,-122.88,25.1128], Tmin=(100,'K'), Tmax=(864.881,'K')), NASAPolynomial(coeffs=[8.05346,0.0272612,-1.25634e-05,2.32479e-09,-1.56266e-13,-812.978,-5.07732], Tmin=(864.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.88283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH]C(C=O)C([CH2])O(13872)',
    structure = SMILES('[CH]C(C=O)C([CH2])O'),
    E0 = (230.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708353,0.0757665,-9.10177e-05,6.04407e-08,-1.63208e-11,27861.2,29.564], Tmin=(100,'K'), Tmax=(897.984,'K')), NASAPolynomial(coeffs=[11.2109,0.028983,-1.28691e-05,2.4221e-09,-1.68083e-13,25975,-19.9745], Tmin=(897.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(O)C([CH2])C=O(13873)',
    structure = SMILES('[CH]C(O)C([CH2])C=O'),
    E0 = (229.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518879,0.0812959,-0.000109701,8.2865e-08,-2.53722e-11,27736.8,29.9581], Tmin=(100,'K'), Tmax=(796.543,'K')), NASAPolynomial(coeffs=[10.6599,0.0303696,-1.37981e-05,2.59698e-09,-1.79018e-13,26121.3,-16.6595], Tmin=(796.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCJ2_triplet)"""),
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
    E0 = (-7.01478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (82.4797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (71.5536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (114.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (89.0488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (63.1984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (75.2828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-7.01478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (193.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (119.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (95.6648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (110.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (110.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (150.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (110.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (69.0098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (76.6652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (387.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (442.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (225.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (333.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (279.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (381.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (363.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (179.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (30.8535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (56.3854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (56.3854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (311.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (376.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (275.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (391.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (346.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (705.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (152.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (115.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (151.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1.26954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (255.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (69.7789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (370.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (413.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (442.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (441.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['C=CO(576)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(O)C1CC1[O](13855)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(89.4945,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 85.9 to 89.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C1C([O])CC1O(13856)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(177207,'s^-1'), n=1.88643, Ea=(78.5684,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 73.3 to 78.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C(O)C(=C)C=O(13857)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH2]C(C=O)C(=C)O(13858)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(578)', 'C=CC=O(5269)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.76206,'m^3/(mol*s)'), n=1.97634, Ea=(9.22116,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=O(373)', '[CH2]C(O)C=C(737)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CO_pri_rad]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CO(576)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.308008,'m^3/(mol*s)'), n=2.06448, Ea=(69.3357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 67.7 to 69.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(D)(132)', '[CH2]C(C=C)C=O(12631)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-CsH_Cds-HH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)C(C)=C[O](13859)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(C=O)[C](C)O(13684)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2][C](O)C(C)C=O(13860)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(O)C(C)[C]=O(13861)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=O)C(C)[O](12750)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(=C[O])C(C)O(13687)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C([C]=O)C(C)O(13688)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C([O])C(C)C=O(13689)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OH(D)(132)', '[CH2][CH]C([CH2])C=O(12931)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2]C([O])C([CH2])C=O(12644)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]O(578)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=O(373)', '[CH2][CH]C([CH2])O(729)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2]C(=C[O])C([CH2])O(13862)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2][C](O)C([CH2])C=O(13863)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH2]C(O)C([CH2])[C]=O(13864)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(O)C1[CH]OC1(13865)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C1[CH]OCC1O(13866)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['C=C(C=O)C(C)O(13694)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['C=C(O)C(C)C=O(13867)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C](C[O])C([CH2])O(9867)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C](O)C([CH2])C[O](13868)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C](O)C([CH2])[CH]O(13869)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([O])C([CH2])C[O](13698)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])C([CH2])[CH]O(13699)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[C-]#[O+](374)', '[CH2]CC([CH2])O(3672)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(O)CC=C[O](12790)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(O)[CH]CC=O(13870)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['[CH2]C(C=O)C[CH]O(12788)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(O)C([CH2])C=O(12792)'],
    products = ['O=CC1CCC1O(12794)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C=COC([CH2])O(12791)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(=CO)C([CH2])O(13871)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(T)(28)', '[CH2]C(O)C=C[O](1110)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CH2(T)(28)', '[CH2]C([CH]O)C=O(13541)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(8)', '[CH]C(C=O)C([CH2])O(13872)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['H(8)', '[CH]C(O)C([CH2])C=O(13873)'],
    products = ['[CH2]C(O)C([CH2])C=O(12792)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '3553',
    isomers = [
        '[CH2]C(O)C([CH2])C=O(12792)',
    ],
    reactants = [
        ('C=CO(576)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3553',
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


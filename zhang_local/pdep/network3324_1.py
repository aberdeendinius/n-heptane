species(
    label = '[CH2]C(C=O)C[CH]C(12590)',
    structure = SMILES('[CH2]C(C=O)C[CH]C'),
    E0 = (124.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,1619.74],'cm^-1')),
        HinderedRotor(inertia=(0.177767,'amu*angstrom^2'), symmetry=1, barrier=(4.08721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17787,'amu*angstrom^2'), symmetry=1, barrier=(4.08958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177882,'amu*angstrom^2'), symmetry=1, barrier=(4.08986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177849,'amu*angstrom^2'), symmetry=1, barrier=(4.0891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1779,'amu*angstrom^2'), symmetry=1, barrier=(4.09026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80362,0.07634,-9.20439e-05,7.62615e-08,-2.70757e-11,15087.1,30.5484], Tmin=(100,'K'), Tmax=(792.621,'K')), NASAPolynomial(coeffs=[5.05372,0.0469085,-2.12382e-05,4.00038e-09,-2.75984e-13,14664.1,12.6137], Tmin=(792.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C=CC(42)',
    structure = SMILES('C=CC'),
    E0 = (6.12372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.597443,'amu*angstrom^2'), symmetry=1, barrier=(13.7364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30977,0.00827491,3.37717e-05,-4.3931e-08,1.58773e-11,767.476,9.64349], Tmin=(100,'K'), Tmax=(988,'K')), NASAPolynomial(coeffs=[5.41204,0.0172866,-6.51359e-06,1.20323e-09,-8.55924e-14,-503.177,-4.80153], Tmin=(988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.12372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = 'C[CH]CC1CC1[O](13101)',
    structure = SMILES('C[CH]CC1CC1[O]'),
    E0 = (214.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21555,0.0489752,5.40429e-06,-3.70914e-08,1.64036e-11,25852.5,27.57], Tmin=(100,'K'), Tmax=(1009.28,'K')), NASAPolynomial(coeffs=[12.2577,0.0331865,-1.27053e-05,2.33254e-09,-1.64119e-13,22198.8,-32.8626], Tmin=(1009.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1CC(C)C1[O](13075)',
    structure = SMILES('[CH2]C1CC(C)C1[O]'),
    E0 = (211.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36454,0.0385688,4.9909e-05,-9.29649e-08,3.87478e-11,25514.3,25.1216], Tmin=(100,'K'), Tmax=(952.023,'K')), NASAPolynomial(coeffs=[15.5792,0.0273074,-8.70459e-06,1.55007e-09,-1.1286e-13,20611.5,-54.2913], Tmin=(952.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(CC(C)OJ)"""),
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
    label = 'C=C(C=O)C[CH]C(13102)',
    structure = SMILES('C=C(C=O)C[CH]C'),
    E0 = (23.3384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.200633,'amu*angstrom^2'), symmetry=1, barrier=(4.61294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568594,'amu*angstrom^2'), symmetry=1, barrier=(13.0731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00642517,'amu*angstrom^2'), symmetry=1, barrier=(4.6056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00640689,'amu*angstrom^2'), symmetry=1, barrier=(4.60149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09483,0.0530058,-2.81709e-05,6.49854e-09,-5.76391e-13,2863.1,24.4584], Tmin=(100,'K'), Tmax=(2488.74,'K')), NASAPolynomial(coeffs=[18.8752,0.0260357,-1.19156e-05,2.1442e-09,-1.38986e-13,-5489.32,-71.7973], Tmin=(2488.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.3384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C=O)C=CC(13103)',
    structure = SMILES('[CH2]C(C=O)C=CC'),
    E0 = (44.4884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,204.531],'cm^-1')),
        HinderedRotor(inertia=(0.280332,'amu*angstrom^2'), symmetry=1, barrier=(8.29129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280229,'amu*angstrom^2'), symmetry=1, barrier=(8.29429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279164,'amu*angstrom^2'), symmetry=1, barrier=(8.29307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280259,'amu*angstrom^2'), symmetry=1, barrier=(8.29228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0337,0.0704402,-7.44412e-05,5.16818e-08,-1.59629e-11,5452.81,26.3682], Tmin=(100,'K'), Tmax=(763.275,'K')), NASAPolynomial(coeffs=[6.36001,0.0425244,-1.95747e-05,3.75475e-09,-2.63433e-13,4639.81,2.11119], Tmin=(763.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.4884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)CC=C(12798)',
    structure = SMILES('[CH2]C(C=O)CC=C'),
    E0 = (57.8685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,180,651.862],'cm^-1')),
        HinderedRotor(inertia=(0.143769,'amu*angstrom^2'), symmetry=1, barrier=(3.30553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48941,'amu*angstrom^2'), symmetry=1, barrier=(11.2525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489405,'amu*angstrom^2'), symmetry=1, barrier=(11.2524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373169,'amu*angstrom^2'), symmetry=1, barrier=(11.2524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12491,0.0657862,-5.35464e-05,2.45267e-08,-4.85209e-12,7061.29,27.1358], Tmin=(100,'K'), Tmax=(1155.11,'K')), NASAPolynomial(coeffs=[9.37644,0.0372122,-1.64409e-05,3.11145e-09,-2.17186e-13,5155.01,-13.8629], Tmin=(1155.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.8685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][CH]C(44)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (279.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.00418548,'amu*angstrom^2'), symmetry=1, barrier=(6.91848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00418537,'amu*angstrom^2'), symmetry=1, barrier=(6.91838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25505,0.0137285,1.00536e-05,-1.43788e-08,4.3875e-12,33590.4,14.1736], Tmin=(100,'K'), Tmax=(1201.86,'K')), NASAPolynomial(coeffs=[3.74312,0.0203097,-8.40105e-06,1.5386e-09,-1.05137e-13,32880.4,9.26373], Tmin=(1201.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJC)"""),
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
    label = 'C=CC[CH]C(122)',
    structure = SMILES('C=CC[CH]C'),
    E0 = (154.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2205.43,2205.43],'cm^-1')),
        HinderedRotor(inertia=(0.178503,'amu*angstrom^2'), symmetry=1, barrier=(9.90382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.50649,'amu*angstrom^2'), symmetry=1, barrier=(28.1001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178511,'amu*angstrom^2'), symmetry=1, barrier=(9.90381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09235,0.03644,-5.91514e-06,-8.60019e-09,3.57549e-12,18612.4,21.9142], Tmin=(100,'K'), Tmax=(1213.84,'K')), NASAPolynomial(coeffs=[6.96893,0.0310549,-1.24641e-05,2.24839e-09,-1.52423e-13,16641.4,-5.79986], Tmin=(1213.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC)"""),
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
    label = 'C[CH]CC(C)=C[O](13104)',
    structure = SMILES('C[CH]CC(C)=C[O]'),
    E0 = (47.7592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,596.689,1448.85],'cm^-1')),
        HinderedRotor(inertia=(0.716177,'amu*angstrom^2'), symmetry=1, barrier=(16.4663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717785,'amu*angstrom^2'), symmetry=1, barrier=(16.5033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00333029,'amu*angstrom^2'), symmetry=1, barrier=(16.553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716986,'amu*angstrom^2'), symmetry=1, barrier=(16.4849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620595,0.0653816,-3.86815e-05,5.88304e-09,1.98237e-12,5873.39,29.2153], Tmin=(100,'K'), Tmax=(1070.72,'K')), NASAPolynomial(coeffs=[13.5936,0.0316474,-1.20589e-05,2.15573e-09,-1.47375e-13,2250.89,-38.2018], Tmin=(1070.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.7592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]CC)C=O(13011)',
    structure = SMILES('[CH2]C([CH]CC)C=O'),
    E0 = (129.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,260.462,3759.83],'cm^-1')),
        HinderedRotor(inertia=(0.00248501,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283218,'amu*angstrom^2'), symmetry=1, barrier=(13.6341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600632,'amu*angstrom^2'), symmetry=1, barrier=(28.9153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121563,'amu*angstrom^2'), symmetry=1, barrier=(5.85208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01106,0.0683705,-5.45666e-05,2.47172e-08,-4.86582e-12,15739.3,29.5665], Tmin=(100,'K'), Tmax=(1155.75,'K')), NASAPolynomial(coeffs=[9.28521,0.0397335,-1.73993e-05,3.27777e-09,-2.28175e-13,13826.8,-11.549], Tmin=(1155.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]CCC([CH2])C=O(13105)',
    structure = SMILES('[CH2]CCC([CH2])C=O'),
    E0 = (135.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,237.582,237.582],'cm^-1')),
        HinderedRotor(inertia=(0.00298656,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239345,'amu*angstrom^2'), symmetry=1, barrier=(9.58686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23935,'amu*angstrom^2'), symmetry=1, barrier=(9.58685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915856,'amu*angstrom^2'), symmetry=1, barrier=(36.6838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46237e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.860591,0.0740324,-7.17913e-05,4.30783e-08,-1.13773e-11,16385.5,29.4113], Tmin=(100,'K'), Tmax=(884.797,'K')), NASAPolynomial(coeffs=[7.67248,0.0432369,-1.95833e-05,3.741e-09,-2.62416e-13,15180.1,-2.61831], Tmin=(884.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C(C)C=O(13106)',
    structure = SMILES('C[CH][CH]C(C)C=O'),
    E0 = (113.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62681,0.057195,-3.41913e-05,1.10962e-08,-1.64697e-12,13782.6,29.7276], Tmin=(100,'K'), Tmax=(1350.48,'K')), NASAPolynomial(coeffs=[6.30404,0.0433414,-1.88039e-05,3.50016e-09,-2.40797e-13,12519.3,5.75741], Tmin=(1350.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCC=O)"""),
)

species(
    label = 'C[CH]CC(C)[C]=O(13107)',
    structure = SMILES('C[CH]CC(C)[C]=O'),
    E0 = (72.7105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38669,0.0619604,-4.6508e-05,2.25137e-08,-5.23652e-12,8835.34,29.2748], Tmin=(100,'K'), Tmax=(940.919,'K')), NASAPolynomial(coeffs=[5.28279,0.0453974,-2.01034e-05,3.80531e-09,-2.65738e-13,8102.16,10.7157], Tmin=(940.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.7105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(=C[O])CCC(13108)',
    structure = SMILES('[CH2]C(=C[O])CCC'),
    E0 = (4.8121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46178,0.062235,-9.26008e-06,-3.77034e-08,2.03724e-11,720.19,27.3285], Tmin=(100,'K'), Tmax=(965.273,'K')), NASAPolynomial(coeffs=[18.5549,0.0244217,-8.2495e-06,1.48364e-09,-1.06794e-13,-4504.11,-68.289], Tmin=(965.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.8121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)CCC(13109)',
    structure = SMILES('[CH2]C([C]=O)CCC'),
    E0 = (88.7749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88994,0.0719683,-6.37144e-05,3.29886e-08,-7.40869e-12,10786.2,28.6687], Tmin=(100,'K'), Tmax=(1031.65,'K')), NASAPolynomial(coeffs=[9.01312,0.0404723,-1.79197e-05,3.39533e-09,-2.37328e-13,9110.12,-10.774], Tmin=(1031.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.7749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2][CH]CC(C)C=O(13110)',
    structure = SMILES('[CH2][CH]CC(C)C=O'),
    E0 = (119.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99312,0.0532087,6.64717e-06,-1.05674e-07,9.81304e-11,14407.5,27.9078], Tmin=(100,'K'), Tmax=(452.084,'K')), NASAPolynomial(coeffs=[4.18323,0.0477687,-2.15487e-05,4.10108e-09,-2.86801e-13,14067,17.5052], Tmin=(452.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C(127)',
    structure = SMILES('[CH2][CH]C[CH]C'),
    E0 = (426.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1956.05,1956.68],'cm^-1')),
        HinderedRotor(inertia=(0.00225683,'amu*angstrom^2'), symmetry=1, barrier=(6.12968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0240973,'amu*angstrom^2'), symmetry=1, barrier=(65.5182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120346,'amu*angstrom^2'), symmetry=1, barrier=(6.13962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120605,'amu*angstrom^2'), symmetry=1, barrier=(6.13892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95832,0.0347404,-1.13952e-05,5.01973e-10,2.04135e-13,51268.6,22.8061], Tmin=(100,'K'), Tmax=(2290.85,'K')), NASAPolynomial(coeffs=[16.749,0.0200389,-7.90959e-06,1.27462e-09,-7.51997e-14,42489.4,-60.5284], Tmin=(2290.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(=C[O])C[CH]C(13111)',
    structure = SMILES('[CH2]C(=C[O])C[CH]C'),
    E0 = (199.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,759.51,3303.02],'cm^-1')),
        HinderedRotor(inertia=(0.0487233,'amu*angstrom^2'), symmetry=1, barrier=(19.9887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868979,'amu*angstrom^2'), symmetry=1, barrier=(19.9795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0491077,'amu*angstrom^2'), symmetry=1, barrier=(20.0067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.58485,'amu*angstrom^2'), symmetry=1, barrier=(82.4227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683082,0.0621781,-2.78066e-05,-1.02323e-08,9.04363e-12,24094.2,29.2697], Tmin=(100,'K'), Tmax=(987.369,'K')), NASAPolynomial(coeffs=[15.3141,0.0265431,-9.58071e-06,1.70815e-09,-1.1886e-13,20052.7,-46.9656], Tmin=(987.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(RCCJC) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH][CH]C)C=O(13112)',
    structure = SMILES('[CH2]C([CH][CH]C)C=O'),
    E0 = (324.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,2483.21,2483.21],'cm^-1')),
        HinderedRotor(inertia=(0.00364868,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18199,'amu*angstrom^2'), symmetry=1, barrier=(5.9664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182014,'amu*angstrom^2'), symmetry=1, barrier=(5.96644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52605,'amu*angstrom^2'), symmetry=1, barrier=(50.0306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00364878,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81186,0.0579324,-1.24455e-05,-8.63356e-08,9.11729e-11,39089.1,29.6164], Tmin=(100,'K'), Tmax=(465.494,'K')), NASAPolynomial(coeffs=[5.61613,0.0426363,-1.92055e-05,3.619e-09,-2.50217e-13,38546.5,12.1478], Tmin=(465.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC([CH2])C=O(12804)',
    structure = SMILES('[CH2][CH]CC([CH2])C=O'),
    E0 = (329.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1839.52],'cm^-1')),
        HinderedRotor(inertia=(0.0901349,'amu*angstrom^2'), symmetry=1, barrier=(2.07238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.090084,'amu*angstrom^2'), symmetry=1, barrier=(2.07121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900575,'amu*angstrom^2'), symmetry=1, barrier=(2.0706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0899999,'amu*angstrom^2'), symmetry=1, barrier=(2.06927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900564,'amu*angstrom^2'), symmetry=1, barrier=(2.07057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802218,0.0775531,-0.000104266,9.08923e-08,-3.2531e-11,39771.4,32.3381], Tmin=(100,'K'), Tmax=(820.937,'K')), NASAPolynomial(coeffs=[5.00864,0.0443281,-2.02997e-05,3.81667e-09,-2.61908e-13,39509.7,15.487], Tmin=(820.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([C]=O)C[CH]C(13113)',
    structure = SMILES('[CH2]C([C]=O)C[CH]C'),
    E0 = (283.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,235.22,1652.09],'cm^-1')),
        HinderedRotor(inertia=(0.131309,'amu*angstrom^2'), symmetry=1, barrier=(5.1554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131308,'amu*angstrom^2'), symmetry=1, barrier=(5.15539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131308,'amu*angstrom^2'), symmetry=1, barrier=(5.1554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131309,'amu*angstrom^2'), symmetry=1, barrier=(5.1554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00304696,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.841883,0.0752333,-9.46957e-05,7.80988e-08,-2.71133e-11,34171.8,31.5676], Tmin=(100,'K'), Tmax=(801.887,'K')), NASAPolynomial(coeffs=[5.93796,0.0423145,-1.90919e-05,3.58272e-09,-2.46347e-13,33595.6,9.61032], Tmin=(801.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)CJ=O) + radical(RCCJC) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C[CH]CC1[CH]OC1(13114)',
    structure = SMILES('C[CH]CC1[CH]OC1'),
    E0 = (198.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2076,0.0530345,-2.26449e-06,-3.81995e-08,2.14268e-11,23946.6,26.1363], Tmin=(100,'K'), Tmax=(861.948,'K')), NASAPolynomial(coeffs=[11.6001,0.0306297,-8.21302e-06,1.1586e-09,-6.96079e-14,21195.8,-28.022], Tmin=(861.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCsJOCs) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1[CH]OC(C)C1(13028)',
    structure = SMILES('[CH2]C1[CH]OC(C)C1'),
    E0 = (117.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30579,0.0413235,4.81729e-05,-1.02728e-07,4.7599e-11,14273.9,21.5918], Tmin=(100,'K'), Tmax=(878.002,'K')), NASAPolynomial(coeffs=[16.628,0.0211805,-2.2582e-06,-1.35042e-11,9.16876e-15,9669.11,-61.2365], Tmin=(878.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C=O)CCC(13115)',
    structure = SMILES('C=C(C=O)CCC'),
    E0 = (-171.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02455,0.0620048,-3.69087e-05,1.01216e-08,-1.09402e-12,-20470,25.6497], Tmin=(100,'K'), Tmax=(2095.6,'K')), NASAPolynomial(coeffs=[19.4416,0.0268504,-1.17454e-05,2.11635e-09,-1.39001e-13,-28188.9,-76.8276], Tmin=(2095.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC=CC(C)C=O(13116)',
    structure = SMILES('CC=CC(C)C=O'),
    E0 = (-166.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48269,0.0589436,-3.5055e-05,1.05093e-08,-1.3413e-12,-19881.3,24.37], Tmin=(100,'K'), Tmax=(1637.02,'K')), NASAPolynomial(coeffs=[9.59879,0.0391121,-1.68833e-05,3.10895e-09,-2.11144e-13,-22538.6,-18.7858], Tmin=(1637.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C=CCC(C)C=O(13117)',
    structure = SMILES('C=CCC(C)C=O'),
    E0 = (-152.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858167,0.0620306,-3.84701e-05,1.18944e-08,-1.49581e-12,-18239.8,27.7582], Tmin=(100,'K'), Tmax=(1803.09,'K')), NASAPolynomial(coeffs=[15.2148,0.0301816,-1.19746e-05,2.098e-09,-1.37522e-13,-23417,-49.9673], Tmin=(1803.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](C[O])C[CH]C(9432)',
    structure = SMILES('[CH2][C](C[O])C[CH]C'),
    E0 = (420.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,185.642,666.472,1533.62,3697.75],'cm^-1')),
        HinderedRotor(inertia=(0.0286572,'amu*angstrom^2'), symmetry=1, barrier=(1.32246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286572,'amu*angstrom^2'), symmetry=1, barrier=(1.32246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286572,'amu*angstrom^2'), symmetry=1, barrier=(1.32246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286572,'amu*angstrom^2'), symmetry=1, barrier=(1.32246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286572,'amu*angstrom^2'), symmetry=1, barrier=(1.32246,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24449,0.066705,-7.26595e-05,6.2552e-08,-2.3503e-11,50658.3,31.981], Tmin=(100,'K'), Tmax=(816.392,'K')), NASAPolynomial(coeffs=[2.05821,0.0510825,-2.25769e-05,4.19665e-09,-2.87093e-13,50913.2,30.5951], Tmin=(816.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCOJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH][CH]C)C[O](13118)',
    structure = SMILES('[CH2]C([CH][CH]C)C[O]'),
    E0 = (462.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,199.427,566.296,3564.95,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0147367,'amu*angstrom^2'), symmetry=1, barrier=(7.41661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147367,'amu*angstrom^2'), symmetry=1, barrier=(7.41661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147367,'amu*angstrom^2'), symmetry=1, barrier=(7.41661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147367,'amu*angstrom^2'), symmetry=1, barrier=(7.41661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147367,'amu*angstrom^2'), symmetry=1, barrier=(7.41661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13872,0.0724568,-9.66599e-05,9.23113e-08,-3.53975e-11,55706.4,34.7198], Tmin=(100,'K'), Tmax=(843.053,'K')), NASAPolynomial(coeffs=[0.561536,0.0534859,-2.42792e-05,4.52934e-09,-3.08658e-13,56575.2,41.9814], Tmin=(843.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(RCCJC) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([CH]O)[CH][CH]C(13119)',
    structure = SMILES('[CH2]C([CH]O)[CH][CH]C'),
    E0 = (416.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,227.903,870.288,3513.51],'cm^-1')),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.028516,'amu*angstrom^2'), symmetry=1, barrier=(0.980342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538474,0.083546,-0.000114004,9.71136e-08,-3.32309e-11,50268.9,35.5643], Tmin=(100,'K'), Tmax=(863.09,'K')), NASAPolynomial(coeffs=[5.82624,0.0439492,-1.89607e-05,3.4423e-09,-2.30527e-13,49918.2,14.0885], Tmin=(863.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
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
    label = '[CH2]CC[CH]C(25)',
    structure = SMILES('[CH2]CC[CH]C'),
    E0 = (231.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2520.08,2520.08],'cm^-1')),
        HinderedRotor(inertia=(0.115943,'amu*angstrom^2'), symmetry=1, barrier=(4.79465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115943,'amu*angstrom^2'), symmetry=1, barrier=(4.79465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32381,'amu*angstrom^2'), symmetry=1, barrier=(54.7441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00289278,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3051.59,'J/mol'), sigma=(5.73385,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.65 K, Pc=36.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95927,0.043164,-1.89531e-05,3.2896e-09,-1.25701e-13,27931,23.7177], Tmin=(100,'K'), Tmax=(1936.29,'K')), NASAPolynomial(coeffs=[12.6506,0.0264639,-1.01885e-05,1.70855e-09,-1.07054e-13,22781,-37.5335], Tmin=(1936.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]CCC=C[O](12588)',
    structure = SMILES('C[CH]CCC=C[O]'),
    E0 = (63.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769553,0.0607664,-2.42806e-05,-9.80463e-09,7.66759e-12,7706.48,30.1673], Tmin=(100,'K'), Tmax=(1023.1,'K')), NASAPolynomial(coeffs=[13.6686,0.0313133,-1.18552e-05,2.14014e-09,-1.48348e-13,3969.12,-37.7236], Tmin=(1023.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]C[CH]CC=O(13120)',
    structure = SMILES('C[CH]C[CH]CC=O'),
    E0 = (119.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22834,0.047793,2.49037e-05,-1.37136e-07,1.22105e-10,14449.2,27.9365], Tmin=(100,'K'), Tmax=(427.851,'K')), NASAPolynomial(coeffs=[3.46105,0.0487361,-2.21137e-05,4.2348e-09,-2.97787e-13,14229.6,21.7023], Tmin=(427.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)C([CH2])C=O(12594)',
    structure = SMILES('[CH2]C(C)C([CH2])C=O'),
    E0 = (129.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0979743,'amu*angstrom^2'), symmetry=1, barrier=(6.95536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09795,'amu*angstrom^2'), symmetry=1, barrier=(6.95581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0979755,'amu*angstrom^2'), symmetry=1, barrier=(6.95551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53315,'amu*angstrom^2'), symmetry=1, barrier=(108.84,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5336,'amu*angstrom^2'), symmetry=1, barrier=(108.845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3619.88,'J/mol'), sigma=(6.36535,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.42 K, Pc=31.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734228,0.0739993,-6.96632e-05,3.92138e-08,-9.42313e-12,15687.9,29.5707], Tmin=(100,'K'), Tmax=(982.602,'K')), NASAPolynomial(coeffs=[9.32769,0.0390164,-1.6259e-05,2.98018e-09,-2.04203e-13,13999.2,-11.7369], Tmin=(982.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(Isobutyl)"""),
)

species(
    label = 'CC1CC(C=O)C1(12595)',
    structure = SMILES('CC1CC(C=O)C1'),
    E0 = (-131.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60327,0.0358881,4.53963e-05,-7.82496e-08,3.05707e-11,-15689.2,23.1854], Tmin=(100,'K'), Tmax=(994.35,'K')), NASAPolynomial(coeffs=[12.4534,0.0333815,-1.28835e-05,2.43351e-09,-1.76086e-13,-19880.8,-39.3256], Tmin=(994.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C=COC[CH]C(12589)',
    structure = SMILES('[CH2]C=COC[CH]C'),
    E0 = (106.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.352258,0.0614689,1.35369e-06,-5.58706e-08,2.87141e-11,12907.9,28.6038], Tmin=(100,'K'), Tmax=(945.842,'K')), NASAPolynomial(coeffs=[21.4131,0.0190232,-5.26873e-06,9.10515e-10,-6.83397e-14,6838.47,-82.8544], Tmin=(945.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(=CO)C[CH]C(13121)',
    structure = SMILES('[CH2]C(=CO)C[CH]C'),
    E0 = (57.7958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,412.373,412.38],'cm^-1')),
        HinderedRotor(inertia=(0.0188627,'amu*angstrom^2'), symmetry=1, barrier=(2.27617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131428,'amu*angstrom^2'), symmetry=1, barrier=(15.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200892,'amu*angstrom^2'), symmetry=1, barrier=(24.2444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131435,'amu*angstrom^2'), symmetry=1, barrier=(15.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131429,'amu*angstrom^2'), symmetry=1, barrier=(15.8601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289674,0.0689416,-3.05276e-05,-1.72578e-08,1.40921e-11,7096.36,29.2714], Tmin=(100,'K'), Tmax=(942.938,'K')), NASAPolynomial(coeffs=[18.446,0.0230594,-7.0733e-06,1.18058e-09,-8.15049e-14,2288,-64.5961], Tmin=(942.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.7958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = 'C[CH]CC=C[O](845)',
    structure = SMILES('C[CH]CC=C[O]'),
    E0 = (86.8142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,905.564,909.361],'cm^-1')),
        HinderedRotor(inertia=(0.568825,'amu*angstrom^2'), symmetry=1, barrier=(13.0784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199356,'amu*angstrom^2'), symmetry=1, barrier=(13.1106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226167,'amu*angstrom^2'), symmetry=1, barrier=(13.0581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52127,0.0447419,-6.09335e-06,-2.3495e-08,1.22038e-11,10539.3,25.2239], Tmin=(100,'K'), Tmax=(982.149,'K')), NASAPolynomial(coeffs=[12.1696,0.0238803,-8.60434e-06,1.54069e-09,-1.07719e-13,7362.16,-31.4821], Tmin=(982.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(32)',
    structure = SMILES('[CH]C'),
    E0 = (351.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,431.535,1804.51],'cm^-1')),
        HinderedRotor(inertia=(0.000906356,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73285,-0.000244773,3.59198e-05,-4.44289e-08,1.65887e-11,42287.5,7.078], Tmin=(100,'K'), Tmax=(940.479,'K')), NASAPolynomial(coeffs=[5.42969,0.0081677,-2.4253e-06,4.22642e-10,-3.09417e-14,41277.1,-4.6789], Tmin=(940.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])C=O(6124)',
    structure = SMILES('[CH2]C([CH2])C=O'),
    E0 = (188.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.235764,'amu*angstrom^2'), symmetry=1, barrier=(5.42068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235141,'amu*angstrom^2'), symmetry=1, barrier=(5.40636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.526268,'amu*angstrom^2'), symmetry=1, barrier=(12.0999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64923,0.0560988,-7.86058e-05,6.57621e-08,-2.23429e-11,22710.5,19.9459], Tmin=(100,'K'), Tmax=(819.334,'K')), NASAPolynomial(coeffs=[6.60175,0.0257593,-1.17819e-05,2.21158e-09,-1.51531e-13,22105.8,-1.6983], Tmin=(819.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CJC(C)C=O)"""),
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
    label = '[CH]CC([CH2])C=O(12942)',
    structure = SMILES('[CH]CC([CH2])C=O'),
    E0 = (402.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180.163,697.062,2610,2611.18],'cm^-1')),
        HinderedRotor(inertia=(0.549407,'amu*angstrom^2'), symmetry=1, barrier=(12.6319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00717067,'amu*angstrom^2'), symmetry=1, barrier=(2.45816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87997,'amu*angstrom^2'), symmetry=1, barrier=(43.2242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549492,'amu*angstrom^2'), symmetry=1, barrier=(12.6339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3289,0.0626831,-7.1191e-05,4.82373e-08,-1.3803e-11,48451.9,24.6975], Tmin=(100,'K'), Tmax=(836.678,'K')), NASAPolynomial(coeffs=[8.00286,0.0307773,-1.39926e-05,2.66328e-09,-1.85943e-13,47335.1,-6.31075], Tmin=(836.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=O)C[C]C(13122)',
    structure = SMILES('[CH2]C(C=O)C[C]C'),
    E0 = (378.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719165,0.0769886,-8.31838e-05,5.36822e-08,-1.4746e-11,45613.1,28.0362], Tmin=(100,'K'), Tmax=(866.373,'K')), NASAPolynomial(coeffs=[8.89371,0.0392482,-1.78436e-05,3.40478e-09,-2.38372e-13,44196.6,-10.229], Tmin=(866.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH]C(C=O)C[CH]C(13123)',
    structure = SMILES('[CH]C(C=O)C[CH]C'),
    E0 = (362.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,201.653,801.651,1068.87,1341.56,1608.84],'cm^-1')),
        HinderedRotor(inertia=(0.152907,'amu*angstrom^2'), symmetry=1, barrier=(3.57398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152907,'amu*angstrom^2'), symmetry=1, barrier=(3.57398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152907,'amu*angstrom^2'), symmetry=1, barrier=(3.57398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152907,'amu*angstrom^2'), symmetry=1, barrier=(3.57398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152907,'amu*angstrom^2'), symmetry=1, barrier=(3.57398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63512,0.0601757,-2.96469e-05,-3.33506e-08,4.23605e-11,43644.1,28.3344], Tmin=(100,'K'), Tmax=(510.633,'K')), NASAPolynomial(coeffs=[5.38654,0.0438111,-1.9827e-05,3.7691e-09,-2.63107e-13,43091.2,11.095], Tmin=(510.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    E0 = (124.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (214.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (211.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (245.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (262.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (272.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (206.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (207.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (124.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (251.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (288.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (283.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (242.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (242.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (268.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (183.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (178.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (369.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (459.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (411.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (536.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (541.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (495.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (310.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (182.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (187.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (187.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (149.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (442.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (525.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (441.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (836.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (284.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (247.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (289.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (132.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (419.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (201.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (502.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (573.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (571.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (590.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (574.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C=CC(42)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]CC1CC1[O](13101)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.93521e+09,'s^-1'), n=0.743095, Ea=(89.4945,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 85.9 to 89.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['[CH2]C1CC(C)C1[O](13075)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(966131,'s^-1'), n=1.86605, Ea=(86.6697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 82.0 to 86.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=C(C=O)C[CH]C(13102)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(72.1434,'m^3/(mol*s)'), n=1.66666, Ea=(10.8177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH2]C(C=O)C=CC(13103)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', '[CH2]C(C=O)CC=C(12798)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]C(44)', 'C=CC=O(5269)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.76206,'m^3/(mol*s)'), n=1.97634, Ea=(9.22116,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-OneDeH_Cds;CJ] + [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=O(373)', 'C=CC[CH]C(122)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC(42)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(28.116,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 25.8 to 28.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]CC(C)=C[O](13104)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH]CC)C=O(13011)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CCC([CH2])C=O(13105)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH][CH]C(C)C=O(13106)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]CC(C)[C]=O(13107)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['[CH2]C(=C[O])CCC(13108)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['[CH2]C([C]=O)CCC(13109)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.164999,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['[CH2][CH]CC(C)C=O(13110)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C(44)', '[CH2]C=C[O](5266)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=O(373)', '[CH2][CH]C[CH]C(127)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2]C(=C[O])C[CH]C(13111)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C([CH][CH]C)C=O(13112)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2][CH]CC([CH2])C=O(12804)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH2]C([C]=O)C[CH]C(13113)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]CC1[CH]OC1(13114)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['[CH2]C1[CH]OC(C)C1(13028)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C=C(C=O)CCC(13115)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['CC=CC(C)C=O(13116)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C=CCC(C)C=O(13117)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C](C[O])C[CH]C(9432)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH][CH]C)C[O](13118)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH]O)[CH][CH]C(13119)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C-]#[O+](374)', '[CH2]CC[CH]C(25)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]CCC=C[O](12588)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['C[CH]C[CH]CC=O(13120)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C)C([CH2])C=O(12594)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=O)C[CH]C(12590)'],
    products = ['CC1CC(C=O)C1(12595)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C=COC[CH]C(12589)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=CO)C[CH]C(13121)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(28)', 'C[CH]CC=C[O](845)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]C(32)', '[CH2]C([CH2])C=O(6124)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH3](11)', '[CH]CC([CH2])C=O(12942)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(8)', '[CH2]C(C=O)C[C]C(13122)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(8)', '[CH]C(C=O)C[CH]C(13123)'],
    products = ['[CH2]C(C=O)C[CH]C(12590)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '3324',
    isomers = [
        '[CH2]C(C=O)C[CH]C(12590)',
    ],
    reactants = [
        ('C=CC(42)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3324',
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


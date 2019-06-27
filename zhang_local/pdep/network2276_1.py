species(
    label = '[O]OC([O])C[CH]O(8162)',
    structure = SMILES('[O]OC([O])C[CH]O'),
    E0 = (-10.8971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180,1871.62],'cm^-1')),
        HinderedRotor(inertia=(0.241447,'amu*angstrom^2'), symmetry=1, barrier=(5.55135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241438,'amu*angstrom^2'), symmetry=1, barrier=(5.55114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241387,'amu*angstrom^2'), symmetry=1, barrier=(5.54997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241435,'amu*angstrom^2'), symmetry=1, barrier=(5.55106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709597,0.0851005,-0.000155547,1.48509e-07,-5.33396e-11,-1204.32,29.6648], Tmin=(100,'K'), Tmax=(865.495,'K')), NASAPolynomial(coeffs=[5.46116,0.0341924,-1.71472e-05,3.25897e-09,-2.21143e-13,-942.593,13.6912], Tmin=(865.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.8971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = '[O]OC=O(5472)',
    structure = SMILES('[O]OC=O'),
    E0 = (-195.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,739.225,739.248,739.254,739.261,739.261],'cm^-1')),
        HinderedRotor(inertia=(0.00030847,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38375,0.00191416,4.59985e-05,-6.61808e-08,2.67351e-11,-23479.8,12.2589], Tmin=(100,'K'), Tmax=(935.456,'K')), NASAPolynomial(coeffs=[10.7708,0.000103189,1.15693e-06,-1.97269e-10,7.46894e-15,-26164.6,-29.8498], Tmin=(935.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-OdOsH) + radical(C(=O)OOJ)"""),
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
    label = '[O]OC([O])CC=O(8069)',
    structure = SMILES('[O]OC([O])CC=O'),
    E0 = (-114.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,2068.27],'cm^-1')),
        HinderedRotor(inertia=(0.25893,'amu*angstrom^2'), symmetry=1, barrier=(5.95332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258464,'amu*angstrom^2'), symmetry=1, barrier=(5.9426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258762,'amu*angstrom^2'), symmetry=1, barrier=(5.94946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42743,0.0684532,-0.000123897,1.25132e-07,-4.76545e-11,-13649.9,26.866], Tmin=(100,'K'), Tmax=(844.392,'K')), NASAPolynomial(coeffs=[2.24964,0.0369954,-1.9051e-05,3.69598e-09,-2.55247e-13,-12806.2,28.8568], Tmin=(844.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])C=CO(8604)',
    structure = SMILES('[O]OC([O])C=CO'),
    E0 = (-113.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.667552,'amu*angstrom^2'), symmetry=1, barrier=(15.3483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667318,'amu*angstrom^2'), symmetry=1, barrier=(15.343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667287,'amu*angstrom^2'), symmetry=1, barrier=(15.3422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968958,0.0675069,-8.47846e-05,5.41261e-08,-1.35722e-11,-13543.9,26.1538], Tmin=(100,'K'), Tmax=(978.594,'K')), NASAPolynomial(coeffs=[13.3275,0.0169903,-7.35077e-06,1.37333e-09,-9.52565e-14,-15962.7,-33.2016], Tmin=(978.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(=O)C[CH]O(8605)',
    structure = SMILES('[O]OC(=O)C[CH]O'),
    E0 = (-256.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14853,0.0643013,-7.77975e-05,4.97005e-08,-1.26781e-11,-30725.8,26.1829], Tmin=(100,'K'), Tmax=(955.852,'K')), NASAPolynomial(coeffs=[11.6353,0.0204169,-8.93042e-06,1.66854e-09,-1.15481e-13,-32730.5,-23.936], Tmin=(955.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(C(=O)OOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]O[O](8201)',
    structure = SMILES('[O][CH]O[O]'),
    E0 = (228.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,1824.75],'cm^-1')),
        HinderedRotor(inertia=(0.270955,'amu*angstrom^2'), symmetry=1, barrier=(6.2298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97073,0.0345271,-8.63149e-05,9.83474e-08,-3.87219e-11,27554.6,13.3771], Tmin=(100,'K'), Tmax=(878.005,'K')), NASAPolynomial(coeffs=[-0.738633,0.0210949,-1.15487e-05,2.23205e-09,-1.51294e-13,29375,37.4477], Tmin=(878.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(ROOJ) + radical(OCJO)"""),
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
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=CC[CH]O(5578)',
    structure = SMILES('O=CC[CH]O'),
    E0 = (-182.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.22553,'amu*angstrom^2'), symmetry=1, barrier=(5.18538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225821,'amu*angstrom^2'), symmetry=1, barrier=(5.19207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225809,'amu*angstrom^2'), symmetry=1, barrier=(5.1918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91437,0.0512399,-7.79859e-05,6.94997e-08,-2.44016e-11,-21922.2,19.1298], Tmin=(100,'K'), Tmax=(848.535,'K')), NASAPolynomial(coeffs=[5.36695,0.024212,-1.11996e-05,2.09419e-09,-1.42355e-13,-22121,5.32112], Tmin=(848.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH)"""),
)

species(
    label = '[O]CCC([O])O[O](2671)',
    structure = SMILES('[O]CCC([O])O[O]'),
    E0 = (34.5102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,2036.23,2037.11],'cm^-1')),
        HinderedRotor(inertia=(0.204609,'amu*angstrom^2'), symmetry=1, barrier=(4.70437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2054,'amu*angstrom^2'), symmetry=1, barrier=(4.72255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204853,'amu*angstrom^2'), symmetry=1, barrier=(4.70998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.299,0.0741491,-0.000138732,1.44464e-07,-5.58615e-11,4233.64,28.8586], Tmin=(100,'K'), Tmax=(851.072,'K')), NASAPolynomial(coeffs=[0.227711,0.0436735,-2.24326e-05,4.33801e-09,-2.986e-13,5702.05,41.4097], Tmin=(851.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.5102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])[CH]CO(8606)',
    structure = SMILES('[O]OC([O])[CH]CO'),
    E0 = (9.22002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1380,1390,370,380,2900,435,180,2717.86,2720.44],'cm^-1')),
        HinderedRotor(inertia=(0.258827,'amu*angstrom^2'), symmetry=1, barrier=(5.95094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258836,'amu*angstrom^2'), symmetry=1, barrier=(5.95115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82414,'amu*angstrom^2'), symmetry=1, barrier=(41.9406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259315,'amu*angstrom^2'), symmetry=1, barrier=(5.96216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.974129,0.0787556,-0.000141629,1.37359e-07,-5.03113e-11,1206.14,31.4704], Tmin=(100,'K'), Tmax=(857.207,'K')), NASAPolynomial(coeffs=[4.2074,0.0358984,-1.80414e-05,3.45036e-09,-2.35714e-13,1672.08,22.321], Tmin=(857.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.22002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCJCOOH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]O[C](O)C[CH]O(8607)',
    structure = SMILES('[O]O[C](O)C[CH]O'),
    E0 = (-31.3559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394545,0.0884478,-0.000153146,1.34772e-07,-4.51792e-11,-3650.12,30.6239], Tmin=(100,'K'), Tmax=(875.285,'K')), NASAPolynomial(coeffs=[9.96658,0.0254452,-1.21725e-05,2.2598e-09,-1.50961e-13,-4588.03,-10.0663], Tmin=(875.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.3559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCsJOH) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[O][C](C[CH]O)OO(2670)',
    structure = SMILES('[O][C](C[CH]O)OO'),
    E0 = (42.3446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,360,370,350,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481576,0.0894499,-0.000161442,1.51807e-07,-5.42771e-11,5207.97,30.4118], Tmin=(100,'K'), Tmax=(851.407,'K')), NASAPolynomial(coeffs=[6.9101,0.0334412,-1.73005e-05,3.34041e-09,-2.2926e-13,5048.67,5.92507], Tmin=(851.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.3446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH) + radical(Cs_P)"""),
)

species(
    label = '[O]O[C]([O])CCO(2677)',
    structure = SMILES('[O]O[C]([O])CCO'),
    E0 = (14.0514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3615,1277.5,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,1804.58,1804.6],'cm^-1')),
        HinderedRotor(inertia=(0.295966,'amu*angstrom^2'), symmetry=1, barrier=(6.80485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295954,'amu*angstrom^2'), symmetry=1, barrier=(6.80456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295968,'amu*angstrom^2'), symmetry=1, barrier=(6.8049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295943,'amu*angstrom^2'), symmetry=1, barrier=(6.8043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99017,0.077414,-0.000135995,1.30205e-07,-4.74323e-11,1787.58,29.7959], Tmin=(100,'K'), Tmax=(853.395,'K')), NASAPolynomial(coeffs=[4.72826,0.0349355,-1.74637e-05,3.34027e-09,-2.28541e-13,2058.37,17.6791], Tmin=(853.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.0514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(O)[CH][CH]O(8608)',
    structure = SMILES('[O]OC(O)[CH][CH]O'),
    E0 = (-36.1873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379356,0.0897799,-0.000158752,1.419e-07,-4.80552e-11,-4231.6,32.2953], Tmin=(100,'K'), Tmax=(877.832,'K')), NASAPolynomial(coeffs=[9.43834,0.026421,-1.27578e-05,2.3717e-09,-1.58285e-13,-4971.32,-5.38297], Tmin=(877.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.1873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCJCOOH) + radical(ROOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C([CH][CH]O)OO(8609)',
    structure = SMILES('[O]C([CH][CH]O)OO'),
    E0 = (37.5132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465162,0.0907961,-0.000167093,1.58984e-07,-5.71662e-11,4626.54,32.0876], Tmin=(100,'K'), Tmax=(855.204,'K')), NASAPolynomial(coeffs=[6.39071,0.0344016,-1.78767e-05,3.45013e-09,-2.36402e-13,4661.8,10.5588], Tmin=(855.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCJCOOH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]CC(O)O[O](8610)',
    structure = SMILES('[O][CH]CC(O)O[O]'),
    E0 = (-10.8971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180,1871.52],'cm^-1')),
        HinderedRotor(inertia=(0.239541,'amu*angstrom^2'), symmetry=1, barrier=(5.50751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240682,'amu*angstrom^2'), symmetry=1, barrier=(5.53375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242917,'amu*angstrom^2'), symmetry=1, barrier=(5.58514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242552,'amu*angstrom^2'), symmetry=1, barrier=(5.57674,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709627,0.0851001,-0.000155545,1.48507e-07,-5.33385e-11,-1204.33,29.6647], Tmin=(100,'K'), Tmax=(865.5,'K')), NASAPolynomial(coeffs=[5.46109,0.0341926,-1.71473e-05,3.25898e-09,-2.21145e-13,-942.566,13.6916], Tmin=(865.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.8971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]CC([O])OO(2672)',
    structure = SMILES('[O][CH]CC([O])OO'),
    E0 = (62.8034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,180,1430.46,1430.54],'cm^-1')),
        HinderedRotor(inertia=(0.181616,'amu*angstrom^2'), symmetry=1, barrier=(4.17572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181548,'amu*angstrom^2'), symmetry=1, barrier=(4.17414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18158,'amu*angstrom^2'), symmetry=1, barrier=(4.17487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25637,'amu*angstrom^2'), symmetry=1, barrier=(51.8784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790097,0.0861891,-0.000164196,1.66091e-07,-6.27195e-11,7654.04,29.4756], Tmin=(100,'K'), Tmax=(849.687,'K')), NASAPolynomial(coeffs=[2.40986,0.0421787,-2.22691e-05,4.33806e-09,-2.99312e-13,8692.23,29.654], Tmin=(849.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.8034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C[CH]O(1051)',
    structure = SMILES('[O][CH]C[CH]O'),
    E0 = (146.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,180,180,1929.54],'cm^-1')),
        HinderedRotor(inertia=(0.197807,'amu*angstrom^2'), symmetry=1, barrier=(4.54797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197798,'amu*angstrom^2'), symmetry=1, barrier=(4.54777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197809,'amu*angstrom^2'), symmetry=1, barrier=(4.54802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53487,0.0654521,-0.000123184,1.19887e-07,-4.31792e-11,17652.5,21.573], Tmin=(100,'K'), Tmax=(883.129,'K')), NASAPolynomial(coeffs=[3.89392,0.0278131,-1.34726e-05,2.50656e-09,-1.67216e-13,18286.9,16.4361], Tmin=(883.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCsJOH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]CC([O])O[O](8061)',
    structure = SMILES('[O][CH]CC([O])O[O]'),
    E0 = (214.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,180,1899.54,1900.29],'cm^-1')),
        HinderedRotor(inertia=(0.222443,'amu*angstrom^2'), symmetry=1, barrier=(5.11441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222983,'amu*angstrom^2'), symmetry=1, barrier=(5.12682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221973,'amu*angstrom^2'), symmetry=1, barrier=(5.1036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05268,0.0825964,-0.000168779,1.74968e-07,-6.61143e-11,25924.5,29.2928], Tmin=(100,'K'), Tmax=(867.121,'K')), NASAPolynomial(coeffs=[0.794527,0.0405663,-2.13067e-05,4.10428e-09,-2.79772e-13,27594.2,39.8709], Tmin=(867.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJOH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])[CH][CH]O(8611)',
    structure = SMILES('[O]OC([O])[CH][CH]O'),
    E0 = (189.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,180,1581.92,1588.38],'cm^-1')),
        HinderedRotor(inertia=(0.279887,'amu*angstrom^2'), symmetry=1, barrier=(6.43515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279604,'amu*angstrom^2'), symmetry=1, barrier=(6.42865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280981,'amu*angstrom^2'), symmetry=1, barrier=(6.46032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280143,'amu*angstrom^2'), symmetry=1, barrier=(6.44105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726306,0.0872264,-0.000171792,1.68078e-07,-6.06948e-11,22897.1,31.9097], Tmin=(100,'K'), Tmax=(875.494,'K')), NASAPolynomial(coeffs=[4.76287,0.0328105,-1.69267e-05,3.21928e-09,-2.17105e-13,23569,20.846], Tmin=(875.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJOH) + radical(CCOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[O]O[C]([O])C[CH]O(8612)',
    structure = SMILES('[O]O[C]([O])C[CH]O'),
    E0 = (194.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,360,370,350,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180,2135.71],'cm^-1')),
        HinderedRotor(inertia=(0.189541,'amu*angstrom^2'), symmetry=1, barrier=(4.35793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189885,'amu*angstrom^2'), symmetry=1, barrier=(4.36582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189323,'amu*angstrom^2'), symmetry=1, barrier=(4.3529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189878,'amu*angstrom^2'), symmetry=1, barrier=(4.36567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.740704,0.0859044,-0.000166227,1.61009e-07,-5.78477e-11,23478.6,30.2411], Tmin=(100,'K'), Tmax=(873.714,'K')), NASAPolynomial(coeffs=[5.2926,0.031832,-1.63398e-05,3.10699e-09,-2.09748e-13,23951.7,16.1544], Tmin=(873.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCsJOH) + radical(ROOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]OC(O)C=CO(8613)',
    structure = SMILES('[O]OC(O)C=CO'),
    E0 = (-339.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244375,0.0746366,-8.83348e-05,5.03436e-08,-1.09055e-11,-40656.4,27.8859], Tmin=(100,'K'), Tmax=(1187.13,'K')), NASAPolynomial(coeffs=[19.1419,0.00865986,-2.06114e-06,2.60778e-10,-1.4447e-14,-44981,-65.8425], Tmin=(1187.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'O=C(C[CH]O)OO(5994)',
    structure = SMILES('O=C(C[CH]O)OO'),
    E0 = (-450.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891949,0.0678404,-7.30232e-05,4.03734e-08,-8.90483e-12,-54094.3,26.3421], Tmin=(100,'K'), Tmax=(1098.56,'K')), NASAPolynomial(coeffs=[13.6785,0.0212838,-9.45496e-06,1.79745e-09,-1.26229e-13,-56903.7,-36.5477], Tmin=(1098.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCsJOH)"""),
)

species(
    label = '[O]OC(=O)CCO(8614)',
    structure = SMILES('[O]OC(=O)CCO'),
    E0 = (-436.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28719,0.0571817,-5.25573e-05,2.54328e-08,-4.98425e-12,-52412.2,26.1297], Tmin=(100,'K'), Tmax=(1220.19,'K')), NASAPolynomial(coeffs=[12.0815,0.0217957,-9.0565e-06,1.66545e-09,-1.14623e-13,-55046.5,-28.0948], Tmin=(1220.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(C(=O)OOJ)"""),
)

species(
    label = '[O]C(C=CO)OO(8615)',
    structure = SMILES('[O]C(C=CO)OO'),
    E0 = (-265.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693921,0.0712532,-8.06888e-05,4.56234e-08,-1.01309e-11,-31814,26.3799], Tmin=(100,'K'), Tmax=(1100.34,'K')), NASAPolynomial(coeffs=[15.3894,0.0178314,-7.86284e-06,1.49973e-09,-1.05826e-13,-35048,-45.9221], Tmin=(1100.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(O)CC=O(8616)',
    structure = SMILES('[O]OC(O)CC=O'),
    E0 = (-339.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10228,0.0707316,-0.000109811,9.74844e-08,-3.43416e-11,-40779.6,27.1746], Tmin=(100,'K'), Tmax=(822.964,'K')), NASAPolynomial(coeffs=[6.855,0.0307303,-1.49562e-05,2.86629e-09,-1.97937e-13,-41318.7,3.0191], Tmin=(822.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = '[O]C(CC=O)OO(2682)',
    structure = SMILES('[O]C(CC=O)OO'),
    E0 = (-266.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17909,0.0718811,-0.00011878,1.15699e-07,-4.41249e-11,-31921,26.9976], Tmin=(100,'K'), Tmax=(812.846,'K')), NASAPolynomial(coeffs=[3.76727,0.0387797,-2.0115e-05,3.95419e-09,-2.76842e-13,-31669,19.1858], Tmin=(812.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + radical(CCOJ)"""),
)

species(
    label = 'O2(S)(5486)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'O[CH]CC1OO1(8617)',
    structure = SMILES('O[CH]CC1OO1'),
    E0 = (-65.9009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803678,0.0563656,-5.59938e-05,2.76303e-08,-5.11824e-12,-7799.02,24.5003], Tmin=(100,'K'), Tmax=(1530,'K')), NASAPolynomial(coeffs=[15.4451,0.00931358,-1.26281e-06,3.44692e-11,3.28937e-15,-11252.4,-49.0068], Tmin=(1530,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.9009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(dioxirane) + radical(CCsJOH)"""),
)

species(
    label = '[O]C1CC(O)O1(8618)',
    structure = SMILES('[O]C1CC(O)O1'),
    E0 = (-271.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37621,0.0296195,6.40622e-06,-2.96482e-08,1.45515e-11,-32538.7,19.789], Tmin=(100,'K'), Tmax=(888.547,'K')), NASAPolynomial(coeffs=[8.27874,0.0199379,-5.76277e-06,8.74976e-10,-5.55182e-14,-34254.4,-11.7418], Tmin=(888.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Oxetane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(O)C([O])O[O](8164)',
    structure = SMILES('[CH2]C(O)C([O])O[O]'),
    E0 = (7.75123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,256.505,1444.44],'cm^-1')),
        HinderedRotor(inertia=(0.240801,'amu*angstrom^2'), symmetry=1, barrier=(5.84524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244344,'amu*angstrom^2'), symmetry=1, barrier=(5.89143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248897,'amu*angstrom^2'), symmetry=1, barrier=(5.88855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163857,'amu*angstrom^2'), symmetry=1, barrier=(5.82883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4658.91,'J/mol'), sigma=(7.47824,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=727.71 K, Pc=25.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.792609,0.0797435,-0.00013446,1.22425e-07,-4.30358e-11,1038.91,30.1066], Tmin=(100,'K'), Tmax=(848.358,'K')), NASAPolynomial(coeffs=[7.19812,0.0309372,-1.52695e-05,2.91169e-09,-1.99181e-13,621.561,4.20288], Tmin=(848.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.75123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = '[O]O(16)',
    structure = SMILES('[O]O'),
    E0 = (-8.19602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1036.72,2034.11,2034.11],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.04595,-0.00173474,1.0377e-05,-1.02207e-08,3.3493e-12,-986.755,4.63579], Tmin=(100,'K'), Tmax=(932.129,'K')), NASAPolynomial(coeffs=[3.21022,0.00367946,-1.27704e-06,2.18051e-10,-1.46343e-14,-910.359,8.18305], Tmin=(932.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.19602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)"""),
)

species(
    label = '[O]C=C[CH]O(8619)',
    structure = SMILES('[O]C=C[CH]O'),
    E0 = (-97.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26548,0.0214751,3.87549e-05,-7.72883e-08,3.4453e-11,-11594.1,19.0133], Tmin=(100,'K'), Tmax=(924.144,'K')), NASAPolynomial(coeffs=[16.9491,0.000981536,2.12395e-06,-4.41897e-10,2.44632e-14,-16146.9,-60.6169], Tmin=(924.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJO)"""),
)

species(
    label = 'O[CH]CC1OOO1(8620)',
    structure = SMILES('O[CH]CC1OOO1'),
    E0 = (-24.6884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30418,0.0454285,-3.31161e-06,-3.30755e-08,1.66892e-11,-2859.97,25.7964], Tmin=(100,'K'), Tmax=(994.524,'K')), NASAPolynomial(coeffs=[16.8342,0.0150888,-6.00043e-06,1.20406e-09,-9.18922e-14,-7537.55,-57.0285], Tmin=(994.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.6884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Cyclobutane) + radical(CCsJOH)"""),
)

species(
    label = '[O]OC1CC(O)O1(8166)',
    structure = SMILES('[O]OC1CC(O)O1'),
    E0 = (-273.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871242,0.0542737,-4.70507e-05,2.18843e-08,-3.9019e-12,-32741.9,25.1614], Tmin=(100,'K'), Tmax=(1601.41,'K')), NASAPolynomial(coeffs=[12.6673,0.0162074,-3.33745e-06,3.32204e-10,-1.36993e-14,-35416.9,-33.8587], Tmin=(1601.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(Oxetane) + radical(ROOJ)"""),
)

species(
    label = '[O]C1CC(O)OO1(8621)',
    structure = SMILES('[O]C1CC(O)OO1'),
    E0 = (-299.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78748,0.0393981,1.12753e-06,-3.38685e-08,1.78374e-11,-35894.4,22.6875], Tmin=(100,'K'), Tmax=(905.949,'K')), NASAPolynomial(coeffs=[12.2223,0.0181405,-4.76226e-06,7.00199e-10,-4.52727e-14,-38803.4,-32.244], Tmin=(905.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-299.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(12dioxolane) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])C[CH]O(8622)',
    structure = SMILES('[O]C([O])C[CH]O'),
    E0 = (-8.70161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,1973.71,1974.81],'cm^-1')),
        HinderedRotor(inertia=(0.154231,'amu*angstrom^2'), symmetry=1, barrier=(3.54608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154807,'amu*angstrom^2'), symmetry=1, barrier=(3.55933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155175,'amu*angstrom^2'), symmetry=1, barrier=(3.56778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.487,0.0689905,-0.00013196,1.36002e-07,-5.17982e-11,-969.274,25.5202], Tmin=(100,'K'), Tmax=(858.835,'K')), NASAPolynomial(coeffs=[1.17363,0.0374323,-1.91747e-05,3.68984e-09,-2.52649e-13,302.242,34.0736], Tmin=(858.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.70161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = '[CH]CC([O])O[O](8222)',
    structure = SMILES('[CH]CC([O])O[O]'),
    E0 = (416.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,180,1497.45,1500.19,1500.26],'cm^-1')),
        HinderedRotor(inertia=(0.194368,'amu*angstrom^2'), symmetry=1, barrier=(4.4689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19308,'amu*angstrom^2'), symmetry=1, barrier=(4.43929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194316,'amu*angstrom^2'), symmetry=1, barrier=(4.46771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62567,0.0599899,-9.96358e-05,9.41361e-08,-3.45507e-11,50213.2,25.128], Tmin=(100,'K'), Tmax=(833.782,'K')), NASAPolynomial(coeffs=[4.75767,0.0288057,-1.44644e-05,2.79222e-09,-1.92971e-13,50252.6,13.9554], Tmin=(833.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])O[O](607)',
    structure = SMILES('[CH2]C([O])O[O]'),
    E0 = (206.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1967.16],'cm^-1')),
        HinderedRotor(inertia=(0.24983,'amu*angstrom^2'), symmetry=1, barrier=(5.74408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249705,'amu*angstrom^2'), symmetry=1, barrier=(5.74121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29491,0.0444229,-7.76836e-05,7.63138e-08,-2.8547e-11,24875.9,21.4728], Tmin=(100,'K'), Tmax=(845.116,'K')), NASAPolynomial(coeffs=[3.56685,0.022769,-1.15016e-05,2.21713e-09,-1.52683e-13,25219.2,18.8536], Tmin=(845.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(CJCOOH)"""),
)

species(
    label = '[CH]O(5471)',
    structure = SMILES('[CH]O'),
    E0 = (205.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,402.686,3356.18],'cm^-1')),
        HinderedRotor(inertia=(0.0105042,'amu*angstrom^2'), symmetry=1, barrier=(23.1306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76003,0.0029575,8.86344e-06,-1.3392e-08,5.33433e-12,24775.7,6.76105], Tmin=(100,'K'), Tmax=(943.117,'K')), NASAPolynomial(coeffs=[5.07489,0.00326005,-9.68482e-07,1.67779e-10,-1.21779e-14,24266.2,-0.891576], Tmin=(943.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[O]O[CH]C[CH]O(8623)',
    structure = SMILES('[O]O[CH]C[CH]O'),
    E0 = (152.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.141051,'amu*angstrom^2'), symmetry=1, barrier=(3.24303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139462,'amu*angstrom^2'), symmetry=1, barrier=(3.20651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142596,'amu*angstrom^2'), symmetry=1, barrier=(3.27855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.477753,'amu*angstrom^2'), symmetry=1, barrier=(10.9845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9859,0.076106,-0.000134411,1.22906e-07,-4.23562e-11,18406,24.8036], Tmin=(100,'K'), Tmax=(885.386,'K')), NASAPolynomial(coeffs=[6.88161,0.0268378,-1.2598e-05,2.31326e-09,-1.53276e-13,18249.1,2.08757], Tmin=(885.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCsJOH)"""),
)

species(
    label = '[O]OC([O])C[C]O(8624)',
    structure = SMILES('[O]OC([O])C[C]O'),
    E0 = (269.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,295.914,296.138,296.225,2066.6],'cm^-1')),
        HinderedRotor(inertia=(0.140472,'amu*angstrom^2'), symmetry=1, barrier=(8.76331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140663,'amu*angstrom^2'), symmetry=1, barrier=(8.76086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140789,'amu*angstrom^2'), symmetry=1, barrier=(8.76021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140651,'amu*angstrom^2'), symmetry=1, barrier=(8.76189,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745577,0.0825532,-0.000150371,1.40638e-07,-4.99046e-11,32555.3,28.7182], Tmin=(100,'K'), Tmax=(854.137,'K')), NASAPolynomial(coeffs=[7.21103,0.0288707,-1.49944e-05,2.89424e-09,-1.98333e-13,32304.6,3.54314], Tmin=(854.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CH2_triplet) + radical(ROOJ)"""),
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
    E0 = (-10.8971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (133.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (105.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1.99837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (93.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-3.02437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-10.8971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (123.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (118.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (147.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (177.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (131.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (64.3836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (121.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (82.6796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (135.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (137.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (426.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (364.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (401.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (406.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (52.5031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (52.5031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (52.5031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (23.6209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (14.0762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (14.0762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-10.8971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (178.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (72.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (166.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (128.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-2.61275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-2.61275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-3.78427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (239.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (450.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (446.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (559.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (481.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['C=CO(576)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[O]OC([O])CC=O(8069)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[O]OC([O])C=CO(8604)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]OC(=O)C[CH]O(8605)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CO(576)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(578)', '[O]OC=O(5472)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O2(2)', 'O=CC[CH]O(5578)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(180.583,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 177.3 to 180.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]CCC([O])O[O](2671)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC([O])[CH]CO(8606)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]O[C](O)C[CH]O(8607)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C](C[CH]O)OO(2670)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]O[C]([O])CCO(2677)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]OC(O)[CH][CH]O(8608)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C([CH][CH]O)OO(8609)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]CC(O)O[O](8610)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(36450.5,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_1;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][CH]CC([O])OO(2672)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.74437e+06,'s^-1'), n=0.972854, Ea=(72.9565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;XH_out] for rate rule [R6HJ_1;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[O][CH]C[CH]O(1051)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[O][CH]CC([O])O[O](8061)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]O(578)', '[O][CH]O[O](8201)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[O]OC([O])[CH][CH]O(8611)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[O]O[C]([O])C[CH]O(8612)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]OC(O)C=CO(8613)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['O=C(C[CH]O)OO(5994)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]OC(=O)CCO(8614)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]C(C=CO)OO(8615)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]OC(O)CC=O(8616)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]C(CC=O)OO(2682)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['O=CC[CH]O(5578)', 'O2(S)(5486)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['O(T)(63)', 'O[CH]CC1OO1(8617)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['O(T)(63)', '[O]C1CC(O)O1(8618)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_sec_rad_intra;OO] for rate rule [R3OO_SS;C_rad/H/NonDeO_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)C([O])O[O](8164)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]O(16)', '[O]C=C[CH]O(8619)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['O[CH]CC1OOO1(8620)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]OC1CC(O)O1(8166)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]OC([O])C[CH]O(8162)'],
    products = ['[O]C1CC(O)OO1(8621)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[O]C([O])C[CH]O(8622)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['OH(D)(132)', '[CH]CC([O])O[O](8222)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([O])O[O](607)', '[CH]O(5471)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O(T)(63)', '[O]O[CH]C[CH]O(8623)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(8)', '[O]OC([O])C[C]O(8624)'],
    products = ['[O]OC([O])C[CH]O(8162)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2276',
    isomers = [
        '[O]OC([O])C[CH]O(8162)',
    ],
    reactants = [
        ('C=CO(576)', '[O]OC=O(5472)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2276',
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


species(
    label = '[CH2][CH]CC[CH]O(853)',
    structure = SMILES('[CH2][CH]CC[CH]O'),
    E0 = (252.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,180.002,351.949,2967.11],'cm^-1')),
        HinderedRotor(inertia=(0.0868877,'amu*angstrom^2'), symmetry=1, barrier=(1.99778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868877,'amu*angstrom^2'), symmetry=1, barrier=(1.99778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868877,'amu*angstrom^2'), symmetry=1, barrier=(1.99778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868877,'amu*angstrom^2'), symmetry=1, barrier=(1.99778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868877,'amu*angstrom^2'), symmetry=1, barrier=(1.99778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04891,0.0709823,-9.482e-05,8.05036e-08,-2.78819e-11,30418.5,29.2159], Tmin=(100,'K'), Tmax=(839.579,'K')), NASAPolynomial(coeffs=[5.59286,0.0383706,-1.69694e-05,3.13389e-09,-2.12613e-13,30041.9,10.3896], Tmin=(839.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = 'C=CCC[CH][O](846)',
    structure = SMILES('C=CCC[CH][O]'),
    E0 = (205.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,325.874,329.126,329.293,1498.65],'cm^-1')),
        HinderedRotor(inertia=(0.0817076,'amu*angstrom^2'), symmetry=1, barrier=(6.16968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0797758,'amu*angstrom^2'), symmetry=1, barrier=(6.15384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0794416,'amu*angstrom^2'), symmetry=1, barrier=(6.18214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31613,0.0457936,7.87204e-06,-1.13111e-07,1.11183e-10,24812.3,21.6845], Tmin=(100,'K'), Tmax=(433.922,'K')), NASAPolynomial(coeffs=[4.64394,0.0388036,-1.7979e-05,3.44718e-09,-2.41256e-13,24474,10.8281], Tmin=(433.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C=CC[CH]O(6458)',
    structure = SMILES('[CH2]C=CC[CH]O'),
    E0 = (119.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.482275,'amu*angstrom^2'), symmetry=1, barrier=(11.0885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680615,'amu*angstrom^2'), symmetry=1, barrier=(19.8253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862193,'amu*angstrom^2'), symmetry=1, barrier=(19.8235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0039146,'amu*angstrom^2'), symmetry=1, barrier=(44.4466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14,0.0602796,-5.214e-05,2.43253e-08,-4.65291e-12,14468.3,24.6829], Tmin=(100,'K'), Tmax=(1239.82,'K')), NASAPolynomial(coeffs=[11.8288,0.0257942,-1.04173e-05,1.89035e-09,-1.2903e-13,11817.9,-29.1821], Tmin=(1239.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]CC=CO(6459)',
    structure = SMILES('[CH2][CH]CC=CO'),
    E0 = (150.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,512.703,512.705],'cm^-1')),
        HinderedRotor(inertia=(0.141788,'amu*angstrom^2'), symmetry=1, barrier=(3.25998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0772763,'amu*angstrom^2'), symmetry=1, barrier=(14.4166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0772858,'amu*angstrom^2'), symmetry=1, barrier=(14.4166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0772957,'amu*angstrom^2'), symmetry=1, barrier=(14.4166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16333,0.0522573,-1.93161e-05,-1.82285e-08,1.28127e-11,18224.2,26.885], Tmin=(100,'K'), Tmax=(935.515,'K')), NASAPolynomial(coeffs=[15.1067,0.0180831,-5.31753e-06,8.67925e-10,-5.9541e-14,14501.9,-45.4049], Tmin=(935.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(RCCJC)"""),
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
    label = '[CH2][CH]CCC[O](852)',
    structure = SMILES('[CH2][CH]CCC[O]'),
    E0 = (297.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1513.37,1513.39,1513.53],'cm^-1')),
        HinderedRotor(inertia=(0.204135,'amu*angstrom^2'), symmetry=1, barrier=(4.69346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203953,'amu*angstrom^2'), symmetry=1, barrier=(4.68928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00288081,'amu*angstrom^2'), symmetry=1, barrier=(4.68043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204178,'amu*angstrom^2'), symmetry=1, barrier=(4.69446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46601,0.0404004,-1.59913e-05,1.98402e-09,-4.37817e-15,35768.9,21.6905], Tmin=(100,'K'), Tmax=(2682.4,'K')), NASAPolynomial(coeffs=[36.4372,0.00376555,-2.51286e-06,3.75858e-10,-1.68211e-14,13572.1,-178.313], Tmin=(2682.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]C[CH]O(6460)',
    structure = SMILES('[CH2]C[CH]C[CH]O'),
    E0 = (252.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,244.623,270.728,3144.33],'cm^-1')),
        HinderedRotor(inertia=(0.0614903,'amu*angstrom^2'), symmetry=1, barrier=(1.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0614903,'amu*angstrom^2'), symmetry=1, barrier=(1.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0614903,'amu*angstrom^2'), symmetry=1, barrier=(1.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0614903,'amu*angstrom^2'), symmetry=1, barrier=(1.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0614903,'amu*angstrom^2'), symmetry=1, barrier=(1.517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08634,0.0714578,-0.00010002,8.89148e-08,-3.17433e-11,30417.3,28.9096], Tmin=(100,'K'), Tmax=(844.671,'K')), NASAPolynomial(coeffs=[4.51083,0.0403032,-1.81681e-05,3.37541e-09,-2.29327e-13,30371.7,16.1208], Tmin=(844.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C[CH]CO(6461)',
    structure = SMILES('[CH2][CH]C[CH]CO'),
    E0 = (271.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,180,1395.16,3549.77],'cm^-1')),
        HinderedRotor(inertia=(0.0333326,'amu*angstrom^2'), symmetry=1, barrier=(1.68898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333326,'amu*angstrom^2'), symmetry=1, barrier=(1.68898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333326,'amu*angstrom^2'), symmetry=1, barrier=(1.68898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333326,'amu*angstrom^2'), symmetry=1, barrier=(1.68898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333326,'amu*angstrom^2'), symmetry=1, barrier=(1.68898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63167,0.0554589,-4.86625e-05,2.92207e-08,-8.16152e-12,32758.2,30.0827], Tmin=(100,'K'), Tmax=(822.616,'K')), NASAPolynomial(coeffs=[5.27163,0.0377592,-1.63876e-05,3.06399e-09,-2.1218e-13,32159.3,13.2328], Tmin=(822.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC[CH][CH]O(6462)',
    structure = SMILES('[CH2]CC[CH][CH]O'),
    E0 = (257.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18257,0.0638329,-6.01213e-05,3.27682e-08,-7.5519e-12,31074.2,28.5035], Tmin=(100,'K'), Tmax=(1024.48,'K')), NASAPolynomial(coeffs=[9.30711,0.0321106,-1.36738e-05,2.5425e-09,-1.75871e-13,29409.5,-10.889], Tmin=(1024.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CCO(6463)',
    structure = SMILES('[CH2][CH][CH]CCO'),
    E0 = (266.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,792.291,2536.58,3733.97],'cm^-1')),
        HinderedRotor(inertia=(0.0286182,'amu*angstrom^2'), symmetry=1, barrier=(5.02145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286182,'amu*angstrom^2'), symmetry=1, barrier=(5.02145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286182,'amu*angstrom^2'), symmetry=1, barrier=(5.02145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286182,'amu*angstrom^2'), symmetry=1, barrier=(5.02145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286182,'amu*angstrom^2'), symmetry=1, barrier=(5.02145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54513,0.0630559,-8.88454e-05,8.61984e-08,-3.28718e-11,32100.8,30.4481], Tmin=(100,'K'), Tmax=(854.671,'K')), NASAPolynomial(coeffs=[0.772468,0.0454075,-2.05507e-05,3.81555e-09,-2.58689e-13,33009.5,38.598], Tmin=(854.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJCC)"""),
)

species(
    label = 'C[CH][CH]C[CH]O(851)',
    structure = SMILES('C[CH][CH]C[CH]O'),
    E0 = (241.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,281.276,1395.15,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0824225,'amu*angstrom^2'), symmetry=1, barrier=(2.22258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0824225,'amu*angstrom^2'), symmetry=1, barrier=(2.22258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0824225,'amu*angstrom^2'), symmetry=1, barrier=(2.22258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0824225,'amu*angstrom^2'), symmetry=1, barrier=(2.22258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0824225,'amu*angstrom^2'), symmetry=1, barrier=(2.22258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26175,0.0707791,-0.000108544,1.04724e-07,-3.88905e-11,29109,29.2282], Tmin=(100,'K'), Tmax=(867.962,'K')), NASAPolynomial(coeffs=[1.50955,0.0446569,-2.02296e-05,3.73312e-09,-2.51196e-13,30007,33.4883], Tmin=(867.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH][CH]O(850)',
    structure = SMILES('C[CH]C[CH][CH]O'),
    E0 = (246.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,279.102,421.651,3036.19],'cm^-1')),
        HinderedRotor(inertia=(0.0734978,'amu*angstrom^2'), symmetry=1, barrier=(1.86218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734978,'amu*angstrom^2'), symmetry=1, barrier=(1.86218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734978,'amu*angstrom^2'), symmetry=1, barrier=(1.86218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734978,'amu*angstrom^2'), symmetry=1, barrier=(1.86218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734978,'amu*angstrom^2'), symmetry=1, barrier=(1.86218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23987,0.0645946,-7.39897e-05,5.61695e-08,-1.83476e-11,29771,29.2434], Tmin=(100,'K'), Tmax=(795.034,'K')), NASAPolynomial(coeffs=[6.24778,0.0365777,-1.58075e-05,2.91854e-09,-1.99343e-13,29063.9,6.79245], Tmin=(795.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CCC[CH][O](848)',
    structure = SMILES('[CH2]CCC[CH][O]'),
    E0 = (283.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,1501.86,1501.9],'cm^-1')),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159195,'amu*angstrom^2'), symmetry=1, barrier=(3.6602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159197,'amu*angstrom^2'), symmetry=1, barrier=(3.66025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159221,'amu*angstrom^2'), symmetry=1, barrier=(3.66081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19727,0.0682175,-8.84193e-05,7.82311e-08,-2.89249e-11,34172.2,26.808], Tmin=(100,'K'), Tmax=(802.529,'K')), NASAPolynomial(coeffs=[4.00256,0.0429162,-1.99724e-05,3.79698e-09,-2.62747e-13,34086.4,16.1622], Tmin=(802.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]CC[CH][O](636)',
    structure = SMILES('C[CH]CC[CH][O]'),
    E0 = (272.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1634.81,1635.01,1635.03],'cm^-1')),
        HinderedRotor(inertia=(0.22758,'amu*angstrom^2'), symmetry=1, barrier=(5.23252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227501,'amu*angstrom^2'), symmetry=1, barrier=(5.2307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227663,'amu*angstrom^2'), symmetry=1, barrier=(5.23443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227522,'amu*angstrom^2'), symmetry=1, barrier=(5.23119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3635.57,'J/mol'), sigma=(6.45581,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.87 K, Pc=30.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35383,0.0677682,-9.77597e-05,9.50667e-08,-3.64621e-11,32864.7,27.1937], Tmin=(100,'K'), Tmax=(841.755,'K')), NASAPolynomial(coeffs=[1.09976,0.0470952,-2.19302e-05,4.12966e-09,-2.82506e-13,33682.6,32.9801], Tmin=(841.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][O](856)',
    structure = SMILES('[CH2][CH]CC[CH][O]'),
    E0 = (477.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,180,1699.8,1701.4,1701.94],'cm^-1')),
        HinderedRotor(inertia=(0.133204,'amu*angstrom^2'), symmetry=1, barrier=(3.06262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132404,'amu*angstrom^2'), symmetry=1, barrier=(3.04423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13148,'amu*angstrom^2'), symmetry=1, barrier=(3.02298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13337,'amu*angstrom^2'), symmetry=1, barrier=(3.06645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3635.57,'J/mol'), sigma=(6.45581,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.87 K, Pc=30.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37832,0.068652,-0.000108722,1.07921e-07,-4.11072e-11,57547.9,28.8922], Tmin=(100,'K'), Tmax=(853.925,'K')), NASAPolynomial(coeffs=[0.96539,0.0446748,-2.10874e-05,3.96918e-09,-2.70396e-13,58563.2,36.3507], Tmin=(853.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]O(6410)',
    structure = SMILES('[CH2][CH][CH]C[CH]O'),
    E0 = (446.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,484.377,2806.56,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0325057,'amu*angstrom^2'), symmetry=1, barrier=(4.70883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325057,'amu*angstrom^2'), symmetry=1, barrier=(4.70883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325057,'amu*angstrom^2'), symmetry=1, barrier=(4.70883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325057,'amu*angstrom^2'), symmetry=1, barrier=(4.70883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325057,'amu*angstrom^2'), symmetry=1, barrier=(4.70883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29139,0.0716024,-0.000119303,1.17346e-07,-4.3461e-11,53792,30.9083], Tmin=(100,'K'), Tmax=(876.806,'K')), NASAPolynomial(coeffs=[1.34208,0.0422944,-1.94209e-05,3.58083e-09,-2.39773e-13,54900.8,37.044], Tmin=(876.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C[CH][CH]O(6409)',
    structure = SMILES('[CH2][CH]C[CH][CH]O'),
    E0 = (451.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,202.333,1873.36,3976.81],'cm^-1')),
        HinderedRotor(inertia=(0.121797,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121797,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121797,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121797,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121797,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23084,0.0659112,-8.6651e-05,7.15117e-08,-2.4189e-11,54455.6,31.0597], Tmin=(100,'K'), Tmax=(837.52,'K')), NASAPolynomial(coeffs=[6.2048,0.0339929,-1.48662e-05,2.73411e-09,-1.85204e-13,53908.8,9.65423], Tmin=(837.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJ) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C=CCCO(6464)',
    structure = SMILES('[CH2]C=CCCO'),
    E0 = (-60.8757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17627,0.054237,-3.01435e-05,3.69312e-09,1.67303e-12,-7213.39,25.0069], Tmin=(100,'K'), Tmax=(1115.22,'K')), NASAPolynomial(coeffs=[11.9559,0.0277474,-1.08889e-05,1.97148e-09,-1.35327e-13,-10374.8,-31.5684], Tmin=(1115.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.8757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CCC=CO(6465)',
    structure = SMILES('[CH2]CCC=CO'),
    E0 = (-43.8483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.926107,0.052499,-1.39432e-06,-4.49379e-08,2.38455e-11,-5149.18,25.001], Tmin=(100,'K'), Tmax=(934.835,'K')), NASAPolynomial(coeffs=[18.4447,0.0157978,-3.8924e-06,6.21333e-10,-4.56505e-14,-10096.3,-67.277], Tmin=(934.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.8483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC=O(6466)',
    structure = SMILES('[CH2]CCCC=O'),
    E0 = (-45.6373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97958,0.0495291,-2.81525e-05,7.88491e-09,-9.37019e-13,-5421.67,22.8927], Tmin=(100,'K'), Tmax=(1703.32,'K')), NASAPolynomial(coeffs=[8.27423,0.034747,-1.51348e-05,2.78984e-09,-1.89199e-13,-7566.02,-10.8276], Tmin=(1703.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.6373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C[CH]O(812)',
    structure = SMILES('[CH2]C([CH2])C[CH]O'),
    E0 = (256.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,276.195,1670.77],'cm^-1')),
        HinderedRotor(inertia=(0.149255,'amu*angstrom^2'), symmetry=1, barrier=(8.06891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149212,'amu*angstrom^2'), symmetry=1, barrier=(8.06877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(6.03862e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23487,'amu*angstrom^2'), symmetry=1, barrier=(66.7863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14935,'amu*angstrom^2'), symmetry=1, barrier=(8.06898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.48678,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=30.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744103,0.0722571,-8.43869e-05,5.74998e-08,-1.57949e-11,31009.2,28.3774], Tmin=(100,'K'), Tmax=(946.422,'K')), NASAPolynomial(coeffs=[10.3085,0.0286393,-1.01934e-05,1.67096e-09,-1.05463e-13,29341.9,-16.4827], Tmin=(946.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC([CH2])O(829)',
    structure = SMILES('[CH2][CH]CC([CH2])O'),
    E0 = (270.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,1165.33,1174.98],'cm^-1')),
        HinderedRotor(inertia=(0.252811,'amu*angstrom^2'), symmetry=1, barrier=(5.81262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258892,'amu*angstrom^2'), symmetry=1, barrier=(5.95244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259807,'amu*angstrom^2'), symmetry=1, barrier=(5.97348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260238,'amu*angstrom^2'), symmetry=1, barrier=(5.98339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00108098,'amu*angstrom^2'), symmetry=1, barrier=(1.05458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.48678,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=30.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24171,0.0641975,-6.8062e-05,4.59783e-08,-1.34321e-11,32657.1,29.2722], Tmin=(100,'K'), Tmax=(814.53,'K')), NASAPolynomial(coeffs=[7.0674,0.0355882,-1.53757e-05,2.85572e-09,-1.96511e-13,31708,2.3618], Tmin=(814.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCO) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1CCC1O(6126)',
    structure = SMILES('[CH2]C1CCC1O'),
    E0 = (13.6632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88137,0.0297487,4.9767e-05,-8.63215e-08,3.56315e-11,1734.95,22.1872], Tmin=(100,'K'), Tmax=(945.144,'K')), NASAPolynomial(coeffs=[13.4828,0.0230825,-6.99694e-06,1.21908e-09,-8.85288e-14,-2353.32,-43.1549], Tmin=(945.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.6632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CCC[CH]O(5624)',
    structure = SMILES('C=CCC[CH]O'),
    E0 = (-19.8316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,322.32,322.469,3961.36],'cm^-1')),
        HinderedRotor(inertia=(0.331115,'amu*angstrom^2'), symmetry=1, barrier=(24.4048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138104,'amu*angstrom^2'), symmetry=1, barrier=(10.1761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241818,'amu*angstrom^2'), symmetry=1, barrier=(17.8052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138016,'amu*angstrom^2'), symmetry=1, barrier=(10.1782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21348,0.0611339,-5.11179e-05,2.38268e-08,-4.67603e-12,-2284.8,24.5777], Tmin=(100,'K'), Tmax=(1188.55,'K')), NASAPolynomial(coeffs=[10.2191,0.0308261,-1.28682e-05,2.37229e-09,-1.63272e-13,-4425.52,-20.4246], Tmin=(1188.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.8316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH)"""),
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
    label = '[CH]CC[CH][CH2](6152)',
    structure = SMILES('[CH]CC[CH][CH2]'),
    E0 = (679.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,220.978,406.057,1393.01,1671.69,3531.85],'cm^-1')),
        HinderedRotor(inertia=(0.0652033,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652033,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652033,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652033,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16391,0.0433952,-2.95389e-05,1.27192e-08,-2.68216e-12,81827.4,23.9734], Tmin=(100,'K'), Tmax=(1002.91,'K')), NASAPolynomial(coeffs=[4.6869,0.0333326,-1.4489e-05,2.71514e-09,-1.8843e-13,81321.4,11.7941], Tmin=(1002.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]O(426)',
    structure = SMILES('[CH2]C[CH]O'),
    E0 = (105.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.100351,'amu*angstrom^2'), symmetry=1, barrier=(2.30727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100273,'amu*angstrom^2'), symmetry=1, barrier=(2.30548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10032,'amu*angstrom^2'), symmetry=1, barrier=(2.30656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27596,0.0395262,-4.11204e-05,2.60331e-08,-6.98322e-12,12712.5,17.592], Tmin=(100,'K'), Tmax=(890.236,'K')), NASAPolynomial(coeffs=[6.59532,0.020118,-8.41765e-06,1.54257e-09,-1.0548e-13,11943.4,-2.74411], Tmin=(890.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH][CH2](721)',
    structure = SMILES('[CH][CH2]'),
    E0 = (556.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1101.59,1101.66],'cm^-1')),
        HinderedRotor(inertia=(0.00420677,'amu*angstrom^2'), symmetry=1, barrier=(3.62356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77493,-0.000462567,3.18167e-05,-4.30783e-08,1.77606e-11,66973.8,8.79001], Tmin=(100,'K'), Tmax=(870.354,'K')), NASAPolynomial(coeffs=[6.06996,0.00332438,5.85464e-07,-2.32999e-10,1.82455e-14,66031.4,-5.08252], Tmin=(870.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCJ)"""),
)

species(
    label = '[CH2][CH]C[CH2](103)',
    structure = SMILES('[CH2][CH]C[CH2]'),
    E0 = (460.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,807.49],'cm^-1')),
        HinderedRotor(inertia=(0.010536,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00642069,'amu*angstrom^2'), symmetry=1, barrier=(45.3249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0169078,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81947,0.0272005,-9.56978e-06,5.91754e-10,1.8328e-13,55442.9,20.1592], Tmin=(100,'K'), Tmax=(1946.61,'K')), NASAPolynomial(coeffs=[9.19054,0.0193114,-7.49953e-06,1.2557e-09,-7.83166e-14,51976.9,-17.353], Tmin=(1946.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]CC[CH]O(6467)',
    structure = SMILES('[CH2][C]CC[CH]O'),
    E0 = (505.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.859336,0.0729834,-9.12933e-05,6.58938e-08,-1.95265e-11,60948.9,27.0739], Tmin=(100,'K'), Tmax=(818.74,'K')), NASAPolynomial(coeffs=[9.55785,0.0304884,-1.34429e-05,2.5065e-09,-1.72323e-13,59524.5,-13.1522], Tmin=(818.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]CC[C]O(6468)',
    structure = SMILES('[CH2][CH]CC[C]O'),
    E0 = (532.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11669,0.0680302,-8.80891e-05,7.04231e-08,-2.34248e-11,64176.8,28.1573], Tmin=(100,'K'), Tmax=(796.44,'K')), NASAPolynomial(coeffs=[7.23782,0.0332363,-1.49285e-05,2.79629e-09,-1.92096e-13,63330.2,0.82618], Tmin=(796.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CH2_triplet) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH]CC[CH]O(6469)',
    structure = SMILES('[CH][CH]CC[CH]O'),
    E0 = (495.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943167,0.0734935,-0.000104527,8.86944e-08,-3.00967e-11,59644.5,28.8129], Tmin=(100,'K'), Tmax=(852.182,'K')), NASAPolynomial(coeffs=[6.77826,0.0344431,-1.5265e-05,2.806e-09,-1.89258e-13,59073.5,4.07973], Tmin=(852.182,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    E0 = (252.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (453.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (339.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (369.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (328.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (349.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (386.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (389.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (380.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (373.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (383.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (396.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (383.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (376.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (366.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (689.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (619.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (658.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (663.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (315.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (315.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (277.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (414.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (429.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (260.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (252.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (713.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (696.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (700.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (717.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (744.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (706.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['[CH2]C=C(87)', 'C=CO(576)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=CCC[CH][O](846)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C=CC[CH]O(6458)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2][CH]CC=CO(6459)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C=C(87)', '[CH2][CH]O(578)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CO(576)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0133609,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]CCC[O](852)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C[CH]C[CH]O(6460)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]C[CH]CO(6461)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC[CH][CH]O(6462)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH]CCO(6463)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['C[CH][CH]C[CH]O(851)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['C[CH]C[CH][CH]O(850)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CCC[CH][O](848)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(72901.1,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH]CC[CH][O](636)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6Hall;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC[CH][O](856)', 'H(8)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]O(578)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH2][CH][CH]C[CH]O(6410)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2][CH]C[CH][CH]O(6409)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['[CH2]C=CCCO(6464)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['[CH2]CCC=CO(6465)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['[CH2]CCCC=O(6466)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH2])C[CH]O(812)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]CC([CH2])O(829)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['[CH2]C1CCC1O(6126)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]CC[CH]O(853)'],
    products = ['C=CCC[CH]O(5624)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OH(D)(132)', '[CH]CC[CH][CH2](6152)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C[CH]O(426)', '[CH][CH2](721)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C[CH2](103)', '[CH]O(5471)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(8)', '[CH2][C]CC[CH]O(6467)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(8)', '[CH2][CH]CC[C]O(6468)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(8)', '[CH][CH]CC[CH]O(6469)'],
    products = ['[CH2][CH]CC[CH]O(853)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '1863',
    isomers = [
        '[CH2][CH]CC[CH]O(853)',
    ],
    reactants = [
        ('[CH2]C=C(87)', 'C=CO(576)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '1863',
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


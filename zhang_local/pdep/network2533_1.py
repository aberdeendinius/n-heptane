species(
    label = '[CH2][CH]CC1([CH2])CO1(9325)',
    structure = SMILES('[CH2][CH]CC1([CH2])CO1'),
    E0 = (417.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,3150,900,1100,180,977.394,977.394,977.394,977.394,977.394,977.394,977.394,977.394,977.394,2360.15],'cm^-1')),
        HinderedRotor(inertia=(0.0738633,'amu*angstrom^2'), symmetry=1, barrier=(1.69826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0738633,'amu*angstrom^2'), symmetry=1, barrier=(1.69826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0738633,'amu*angstrom^2'), symmetry=1, barrier=(1.69826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0738633,'amu*angstrom^2'), symmetry=1, barrier=(1.69826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.221612,0.0813842,-9.42854e-05,6.21941e-08,-1.60021e-11,50344.9,29.0871], Tmin=(100,'K'), Tmax=(1088.93,'K')), NASAPolynomial(coeffs=[11.88,0.0286761,-8.06658e-06,1.07456e-09,-5.66603e-14,48391.8,-25.4609], Tmin=(1088.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(RCCJ) + radical(CJC(C)OC)"""),
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
    label = 'C=C1CO1(9173)',
    structure = SMILES('C=C1CO1'),
    E0 = (24.285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,3150,900,1100,593.567,593.583,593.593,593.632,593.635,593.636,593.645,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3128.07,'J/mol'), sigma=(5.27733,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=488.60 K, Pc=48.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97478,0.00612362,6.24447e-05,-9.43839e-08,3.93209e-11,2973,8.45963], Tmin=(100,'K'), Tmax=(918.038,'K')), NASAPolynomial(coeffs=[14.1102,0.00013146,2.75083e-06,-5.76257e-10,3.43139e-14,-863.606,-54.0705], Tmin=(918.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane)"""),
)

species(
    label = '[CH2][CH]CC(=C)O[CH2](9949)',
    structure = SMILES('[CH2][CH]CC(=C)O[CH2]'),
    E0 = (355.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722819,0.0656633,-5.09849e-05,2.0894e-08,-3.4945e-12,42882.8,32.7895], Tmin=(100,'K'), Tmax=(1408.29,'K')), NASAPolynomial(coeffs=[13.9543,0.0280813,-1.09552e-05,1.94421e-09,-1.30505e-13,39156.1,-35.5747], Tmin=(1408.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=COCJ) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC(=C)C[O](9950)',
    structure = SMILES('[CH2][CH]CC(=C)C[O]'),
    E0 = (392.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,180,1411.66,1412.4,1413.46],'cm^-1')),
        HinderedRotor(inertia=(0.00299583,'amu*angstrom^2'), symmetry=1, barrier=(4.24569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183624,'amu*angstrom^2'), symmetry=1, barrier=(4.22188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184454,'amu*angstrom^2'), symmetry=1, barrier=(4.24096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18465,'amu*angstrom^2'), symmetry=1, barrier=(4.24547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2088,0.0466794,-2.11287e-05,3.85344e-09,-2.49663e-13,47252,24.854], Tmin=(100,'K'), Tmax=(3067.95,'K')), NASAPolynomial(coeffs=[50.8233,-0.00781427,1.80567e-06,-3.24223e-10,2.50876e-14,14466,-264.053], Tmin=(3067.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCOJ) + radical(RCCJC)"""),
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
    label = '[CH2]C1([CH]C=C)CO1(9951)',
    structure = SMILES('[CH2]C1([CH]C=C)CO1'),
    E0 = (262.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,3150,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,1039.94,1039.94,1039.94,1039.94,1039.94,1039.94,1039.94,1039.94,1039.94,2284.68],'cm^-1')),
        HinderedRotor(inertia=(0.078612,'amu*angstrom^2'), symmetry=1, barrier=(1.80744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078612,'amu*angstrom^2'), symmetry=1, barrier=(1.80744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078612,'amu*angstrom^2'), symmetry=1, barrier=(1.80744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.485006,0.0832147,-8.56354e-05,4.55682e-08,-9.1556e-12,31740.7,24.7737], Tmin=(100,'K'), Tmax=(1415.37,'K')), NASAPolynomial(coeffs=[18.1572,0.0177206,-2.65038e-06,8.66165e-11,7.29959e-15,27746.5,-67.1081], Tmin=(1415.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CCJCO) + radical(CJC(C)OC)"""),
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
    label = '[CH2][C]1CO1(9171)',
    structure = SMILES('[CH2][C]1CO1'),
    E0 = (288.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3100,440,815,1455,1000,180,180,1058.34,1058.35,1058.51,1058.53,1058.6],'cm^-1')),
        HinderedRotor(inertia=(0.171032,'amu*angstrom^2'), symmetry=1, barrier=(3.93236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31186,0.0306339,-2.96933e-05,1.70245e-08,-3.60957e-12,34736.5,14.5514], Tmin=(100,'K'), Tmax=(1467.19,'K')), NASAPolynomial(coeffs=[6.04902,0.0115269,-1.04116e-06,-1.37533e-10,2.06889e-14,34599.8,-1.63964], Tmin=(1467.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C[CH]C1([CH2])CO1(9952)',
    structure = SMILES('[CH2]C[CH]C1([CH2])CO1'),
    E0 = (422.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,3025,407.5,1350,352.5,180,1009.57,1009.57,1009.57,1009.57,1009.57,1009.57,1009.57,1009.57,1009.57,2316.75],'cm^-1')),
        HinderedRotor(inertia=(0.0857162,'amu*angstrom^2'), symmetry=1, barrier=(1.97078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0857162,'amu*angstrom^2'), symmetry=1, barrier=(1.97078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0857162,'amu*angstrom^2'), symmetry=1, barrier=(1.97078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0857162,'amu*angstrom^2'), symmetry=1, barrier=(1.97078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.167738,0.0805472,-8.24379e-05,4.54388e-08,-9.59512e-12,51023.2,30.242], Tmin=(100,'K'), Tmax=(1330.25,'K')), NASAPolynomial(coeffs=[15.843,0.0218816,-4.42125e-06,3.93907e-10,-1.21246e-14,47694.5,-48.07], Tmin=(1330.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(RCCJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH][CH]C1(C)CO1(9953)',
    structure = SMILES('[CH2][CH][CH]C1(C)CO1'),
    E0 = (406.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.368526,0.0680779,-5.83455e-05,2.87894e-08,-5.60316e-12,49072.3,31.3061], Tmin=(100,'K'), Tmax=(1431.22,'K')), NASAPolynomial(coeffs=[12.4436,0.0264535,-6.46561e-06,7.78226e-10,-3.8587e-14,46422.6,-28.4602], Tmin=(1431.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2][CH]CC1(C)[CH]O1(9954)',
    structure = SMILES('[CH2][CH]CC1(C)[CH]O1'),
    E0 = (390.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0586942,0.0759056,-7.31725e-05,3.90882e-08,-8.10121e-12,47136.7,30.1203], Tmin=(100,'K'), Tmax=(1345.39,'K')), NASAPolynomial(coeffs=[14.3991,0.0242661,-5.56021e-06,6.10826e-10,-2.7058e-14,44092.9,-40.2901], Tmin=(1345.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(RCCJ) + radical(CCsJO)"""),
)

species(
    label = '[CH2]C1([CH][CH]C)CO1(9425)',
    structure = SMILES('[CH2]C1([CH][CH]C)CO1'),
    E0 = (412.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,3150,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,972.122,972.122,972.122,972.122,972.122,972.122,972.122,972.122,972.122,2361.43],'cm^-1')),
        HinderedRotor(inertia=(0.0711661,'amu*angstrom^2'), symmetry=1, barrier=(1.63625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0711661,'amu*angstrom^2'), symmetry=1, barrier=(1.63625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0711661,'amu*angstrom^2'), symmetry=1, barrier=(1.63625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0711661,'amu*angstrom^2'), symmetry=1, barrier=(1.63625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.243113,0.0770826,-8.11783e-05,4.85179e-08,-1.13073e-11,49704.7,29.7167], Tmin=(100,'K'), Tmax=(1195.47,'K')), NASAPolynomial(coeffs=[12.7138,0.0265352,-6.68674e-06,8.04796e-10,-3.86992e-14,47353.3,-30.0373], Tmin=(1195.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(RCCJC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]CCC1([CH2])[CH]O1(9955)',
    structure = SMILES('[CH2]CCC1([CH2])[CH]O1'),
    E0 = (406.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.504081,0.0886567,-9.81169e-05,5.6672e-08,-1.24256e-11,49088.8,29.1535], Tmin=(100,'K'), Tmax=(1279.87,'K')), NASAPolynomial(coeffs=[17.8304,0.0196302,-3.47596e-06,2.16696e-10,2.36982e-16,45356,-60.072], Tmin=(1279.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CJC(C)OC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1([CH]O1)C[CH]C(9426)',
    structure = SMILES('[CH2]C1([CH]O1)C[CH]C'),
    E0 = (395.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2950,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.121116,0.0855019,-9.78548e-05,6.09266e-08,-1.45894e-11,47771.5,28.7296], Tmin=(100,'K'), Tmax=(1166.79,'K')), NASAPolynomial(coeffs=[14.6899,0.0242871,-5.73771e-06,6.25833e-10,-2.61428e-14,45025.8,-41.9648], Tmin=(1166.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(RCCJC) + radical(CCsJO)"""),
)

species(
    label = '[CH2][CH][CH]C1([CH2])CO1(9956)',
    structure = SMILES('[CH2][CH][CH]C1([CH2])CO1'),
    E0 = (617.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3000,3050,390,425,1340,1360,335,370,180,915.018,915.018,915.018,915.018,915.018,915.018,915.018,915.018,915.018,2410.66],'cm^-1')),
        HinderedRotor(inertia=(0.0468249,'amu*angstrom^2'), symmetry=1, barrier=(1.0766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468249,'amu*angstrom^2'), symmetry=1, barrier=(1.0766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468249,'amu*angstrom^2'), symmetry=1, barrier=(1.0766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468249,'amu*angstrom^2'), symmetry=1, barrier=(1.0766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.352138,0.0769639,-8.86268e-05,5.68511e-08,-1.4065e-11,74384.2,31.1125], Tmin=(100,'K'), Tmax=(1142.16,'K')), NASAPolynomial(coeffs=[12.3726,0.0244806,-6.06023e-06,6.96284e-10,-3.09517e-14,72315.8,-25.5112], Tmin=(1142.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(CJC(C)OC) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC1([CH2])[CH]O1(9957)',
    structure = SMILES('[CH2][CH]CC1([CH2])[CH]O1'),
    E0 = (601.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.021972,0.0854963,-0.000105683,6.97275e-08,-1.75354e-11,72451.5,30.1611], Tmin=(100,'K'), Tmax=(1124.55,'K')), NASAPolynomial(coeffs=[14.3655,0.0222025,-5.09348e-06,5.13072e-10,-1.80395e-14,69981.8,-37.5324], Tmin=(1124.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(RCCJ) + radical(CJC(C)OC) + radical(RCCJC)"""),
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
    label = '[CH2][CH]CC(=C)[O](4399)',
    structure = SMILES('[CH2][CH]CC(=C)[O]'),
    E0 = (282.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,180,467.081,3411.18],'cm^-1')),
        HinderedRotor(inertia=(0.885736,'amu*angstrom^2'), symmetry=1, barrier=(20.3648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000482417,'amu*angstrom^2'), symmetry=1, barrier=(3.9835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104638,'amu*angstrom^2'), symmetry=1, barrier=(84.4336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67971,0.0499013,-4.01055e-05,1.84346e-08,-3.58897e-12,34077.8,26.7309], Tmin=(100,'K'), Tmax=(1197.46,'K')), NASAPolynomial(coeffs=[8.80394,0.0261036,-1.02955e-05,1.8384e-09,-1.24118e-13,32371.6,-8.92326], Tmin=(1197.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'O(S)(1732)',
    structure = SMILES('O'),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,51997.4,2.99252], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,51997.4,2.99252], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][CH]CC([CH2])=C(6303)',
    structure = SMILES('[CH2][CH]CC([CH2])=C'),
    E0 = (471.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,419.284,2014.23],'cm^-1')),
        HinderedRotor(inertia=(0.0636814,'amu*angstrom^2'), symmetry=1, barrier=(7.92145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0634849,'amu*angstrom^2'), symmetry=1, barrier=(7.92449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0171016,'amu*angstrom^2'), symmetry=1, barrier=(30.2503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0489515,'amu*angstrom^2'), symmetry=1, barrier=(86.5653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42509,0.0531145,-3.31274e-05,1.07407e-08,-1.45659e-12,56843.8,26.4331], Tmin=(100,'K'), Tmax=(1641.34,'K')), NASAPolynomial(coeffs=[10.9943,0.029794,-1.1815e-05,2.08414e-09,-1.38063e-13,53702.6,-24.4744], Tmin=(1641.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]C1(C)CO1(9958)',
    structure = SMILES('C=C[CH]C1(C)CO1'),
    E0 = (51.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622525,0.0613309,-9.15137e-06,-4.37682e-08,2.61677e-11,6381.35,21.0604], Tmin=(100,'K'), Tmax=(881.31,'K')), NASAPolynomial(coeffs=[17.8018,0.0207772,-3.81388e-06,3.6892e-10,-1.82039e-14,1900.17,-67.8931], Tmin=(881.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][CH]CC[C]1CO1(9324)',
    structure = SMILES('[CH2][CH]CC[C]1CO1'),
    E0 = (405.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,3150,900,1100,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966438,0.0681907,-7.19537e-05,4.99507e-08,-1.42714e-11,48821.2,30.0137], Tmin=(100,'K'), Tmax=(987.059,'K')), NASAPolynomial(coeffs=[6.76599,0.0369596,-1.27479e-05,2.03003e-09,-1.25005e-13,48052.8,4.01694], Tmin=(987.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(C2CsJO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]1CCO1(9959)',
    structure = SMILES('[CH2][CH]C[C]1CCO1'),
    E0 = (400.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27316,0.0598332,-5.13176e-05,3.12315e-08,-8.3219e-12,48218.3,28.8065], Tmin=(100,'K'), Tmax=(1001.72,'K')), NASAPolynomial(coeffs=[5.80685,0.0380227,-1.31074e-05,2.10765e-09,-1.31514e-13,47496,7.85459], Tmin=(1001.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJC) + radical(C2CsJOCs) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C1([CH2])CO1(9327)',
    structure = SMILES('[CH2]C([CH2])C1([CH2])CO1'),
    E0 = (418.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1380,1390,370,380,2900,435,2750,3150,900,1100,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3609.07,'J/mol'), sigma=(6.55841,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.73 K, Pc=29.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.370629,0.0862578,-9.75577e-05,5.88769e-08,-1.3288e-11,50545.4,29.267], Tmin=(100,'K'), Tmax=(1298.86,'K')), NASAPolynomial(coeffs=[15.2208,0.0210631,-2.42759e-06,-1.3372e-10,3.02046e-14,47944.3,-44.4512], Tmin=(1298.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC2(CO2)C1(9328)',
    structure = SMILES('[CH2]C1CC2(CO2)C1'),
    E0 = (172.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38166,0.0473155,1.10179e-06,-3.02746e-08,1.39267e-11,20822.3,21.4699], Tmin=(100,'K'), Tmax=(995.362,'K')), NASAPolynomial(coeffs=[10.8499,0.0326592,-1.20644e-05,2.15539e-09,-1.48966e-13,17778.6,-29.9859], Tmin=(995.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + polycyclic(s1_3_4_ane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])O[CH2](9960)',
    structure = SMILES('[CH2][CH]C[C]([CH2])O[CH2]'),
    E0 = (656.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,360,370,350,200.451,809.89,1204.92,1618.4],'cm^-1')),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154262,'amu*angstrom^2'), symmetry=1, barrier=(3.5651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0959789,0.105859,-0.000187746,1.77421e-07,-6.28878e-11,79124.2,33.8823], Tmin=(100,'K'), Tmax=(883.293,'K')), NASAPolynomial(coeffs=[5.12643,0.044961,-2.10732e-05,3.88297e-09,-2.58351e-13,79654.7,17.5606], Tmin=(883.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(RCCJC) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(RCCJ) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2][CH]CC([CH2])([CH2])[O](9961)',
    structure = SMILES('[CH2][CH]CC([CH2])([CH2])[O]'),
    E0 = (673.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0514656,0.0975582,-0.000147496,1.2595e-07,-4.24344e-11,81161.5,32.3285], Tmin=(100,'K'), Tmax=(852.257,'K')), NASAPolynomial(coeffs=[9.05467,0.0393747,-1.7908e-05,3.31796e-09,-2.24258e-13,80170.3,-6.85686], Tmin=(852.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CJC(C)2O) + radical(RCCJ) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])C[O](9962)',
    structure = SMILES('[CH2][CH]C[C]([CH2])C[O]'),
    E0 = (625.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,360,370,350,667.488,985.502,3806.13,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0160245,'amu*angstrom^2'), symmetry=1, barrier=(7.39761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160245,'amu*angstrom^2'), symmetry=1, barrier=(7.39761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160245,'amu*angstrom^2'), symmetry=1, barrier=(7.39761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160245,'amu*angstrom^2'), symmetry=1, barrier=(7.39761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160245,'amu*angstrom^2'), symmetry=1, barrier=(7.39761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25486,0.0677675,-8.43021e-05,7.63621e-08,-2.85836e-11,75342.1,33.7293], Tmin=(100,'K'), Tmax=(837.367,'K')), NASAPolynomial(coeffs=[1.97203,0.0485761,-2.16828e-05,4.02374e-09,-2.73933e-13,75774.8,33.6971], Tmin=(837.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ) + radical(CCJ(C)CO) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1(CC=C)CO1(9322)',
    structure = SMILES('[CH2]C1(CC=C)CO1'),
    E0 = (145.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.234544,0.0788709,-7.64837e-05,3.9849e-08,-7.94066e-12,17668.8,26.6764], Tmin=(100,'K'), Tmax=(1419.87,'K')), NASAPolynomial(coeffs=[16.1256,0.0216025,-4.17326e-06,3.52116e-10,-9.99434e-15,14149.8,-54.0187], Tmin=(1419.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CJC(C)OC)"""),
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
    label = '[CH2][CH]C[C]1CO1(6530)',
    structure = SMILES('[CH2][CH]C[C]1CO1'),
    E0 = (428.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,3025,407.5,1350,352.5,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64812,0.0530058,-5.67547e-05,4.01898e-08,-1.14273e-11,51657,25.3206], Tmin=(100,'K'), Tmax=(1038.26,'K')), NASAPolynomial(coeffs=[5.53719,0.0290673,-9.23223e-06,1.36807e-09,-7.91925e-14,51332.1,8.7365], Tmin=(1038.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C1([CH2])CO1(9287)',
    structure = SMILES('[CH2]C1([CH2])CO1'),
    E0 = (275.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867653,0.0644443,-8.34786e-05,5.48812e-08,-1.32592e-11,33291,17.5002], Tmin=(100,'K'), Tmax=(1227.31,'K')), NASAPolynomial(coeffs=[12.5263,0.0113524,-1.42155e-07,-4.076e-10,4.43033e-14,31566,-36.5028], Tmin=(1227.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]CC1([CH2])CO1(9963)',
    structure = SMILES('[CH2][C]CC1([CH2])CO1'),
    E0 = (671.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194685,0.0862309,-0.000101574,6.29736e-08,-1.48763e-11,80885,27.7475], Tmin=(100,'K'), Tmax=(1183.36,'K')), NASAPolynomial(coeffs=[16.0843,0.0203262,-4.2458e-06,3.73471e-10,-1.00006e-14,77793.9,-50.3116], Tmin=(1183.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJ) + radical(CJC(C)OC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1(C[CH][CH2])CO1(9964)',
    structure = SMILES('[CH]C1(C[CH][CH2])CO1'),
    E0 = (655.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401505,0.0731124,-7.54648e-05,4.40409e-08,-1.00436e-11,78931.5,29.6911], Tmin=(100,'K'), Tmax=(1213.9,'K')), NASAPolynomial(coeffs=[12.6098,0.0250444,-6.38026e-06,7.79845e-10,-3.83309e-14,76545.1,-29.1941], Tmin=(1213.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJ2_triplet) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH]CC1([CH2])CO1(9965)',
    structure = SMILES('[CH][CH]CC1([CH2])CO1'),
    E0 = (660.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126762,0.0837807,-0.000103679,7.01552e-08,-1.8224e-11,79570.5,28.6442], Tmin=(100,'K'), Tmax=(1091.03,'K')), NASAPolynomial(coeffs=[12.9102,0.0250098,-6.51171e-06,7.81808e-10,-3.62102e-14,77489.5,-30.8951], Tmin=(1091.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJC) + radical(CCJ2_triplet) + radical(CJC(C)OC)"""),
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
    E0 = (417.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (419.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (479.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (482.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (529.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (481.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (559.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (535.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (534.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (558.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (511.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (472.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (772.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (829.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (812.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (705.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (904.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (480.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (574.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (574.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (576.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (425.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (662.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (679.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (628.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (417.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (844.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (866.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (883.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (866.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (872.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]C=C(87)', 'C=C1CO1(9173)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC(=C)O[CH2](9949)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(64.0574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]CC(=C)C[O](9950)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.33596e+11,'s^-1'), n=-0.0500183, Ea=(86.0977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH2]C1([CH]C=C)CO1(9951)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH2](6136)', 'C=C1CO1(9173)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00359021,'m^3/(mol*s)'), n=2.50446, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CJ] for rate rule [Cds-OsCs_Cds;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C=C(87)', '[CH2][C]1CO1(9171)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C[CH]C1([CH2])CO1(9952)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2][CH][CH]C1(C)CO1(9953)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2][CH]CC1(C)[CH]O1(9954)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;Cs_H_out_H/NonDeO] for rate rule [R3H_SS_23cy3;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]C1([CH][CH]C)CO1(9425)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]CCC1([CH2])[CH]O1(9955)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.22366e+06,'s^-1'), n=1.605, Ea=(93.8055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]C1([CH]O1)C[CH]C(9426)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(126713,'s^-1'), n=1.75034, Ea=(55.0506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]1CO1(9171)', '[CH2][CH][CH2](6136)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', '[CH2][CH][CH]C1([CH2])CO1(9956)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2][CH]CC1([CH2])[CH]O1(9957)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(S)(14)', '[CH2][CH]CC(=C)[O](4399)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.23836e+08,'m^3/(mol*s)'), n=-0.586333, Ea=(3.56505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;multiplebond] for rate rule [carbene;mb_carbonyl]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O(S)(1732)', '[CH2][CH]CC([CH2])=C(6303)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.19311e+07,'m^3/(mol*s)'), n=0, Ea=(0.08368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [o_atom_singlet;mb_db]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['C=C[CH]C1(C)CO1(9958)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2][CH]CC[C]1CO1(9324)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2][CH]C[C]1CCO1(9959)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH2])C1([CH2])CO1(9327)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]C1CC2(CO2)C1(9328)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C[C]([CH2])O[CH2](9960)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]CC([CH2])([CH2])[O](9961)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.18842e+14,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C[C]([CH2])C[O](9962)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    products = ['[CH2]C1(CC=C)CO1(9322)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(T)(28)', '[CH2][CH]C[C]1CO1(6530)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/NDMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH][CH2](721)', '[CH2]C1([CH2])CO1(9287)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(8)', '[CH2][C]CC1([CH2])CO1(9963)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(8)', '[CH]C1(C[CH][CH2])CO1(9964)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(8)', '[CH][CH]CC1([CH2])CO1(9965)'],
    products = ['[CH2][CH]CC1([CH2])CO1(9325)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2533',
    isomers = [
        '[CH2][CH]CC1([CH2])CO1(9325)',
    ],
    reactants = [
        ('[CH2]C=C(87)', 'C=C1CO1(9173)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2533',
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


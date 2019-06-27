species(
    label = '[CH2]CC[CH]O(5497)',
    structure = SMILES('[CH2]CC[CH]O'),
    E0 = (81.4127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,464.832,1913.71],'cm^-1')),
        HinderedRotor(inertia=(0.0116213,'amu*angstrom^2'), symmetry=1, barrier=(1.78982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.443866,'amu*angstrom^2'), symmetry=1, barrier=(10.2054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325503,'amu*angstrom^2'), symmetry=1, barrier=(10.2047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0663452,'amu*angstrom^2'), symmetry=1, barrier=(10.2011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65868,0.0539342,-5.3536e-05,3.21215e-08,-8.24506e-12,9873.9,22.0555], Tmin=(100,'K'), Tmax=(922.554,'K')), NASAPolynomial(coeffs=[7.51718,0.0285328,-1.22349e-05,2.27585e-09,-1.57211e-13,8792.96,-5.7361], Tmin=(922.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.4127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=C(18)',
    structure = SMILES('C=C'),
    E0 = (42.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97974,-0.0075756,5.52973e-05,-6.36222e-08,2.31767e-11,5077.46,4.04623], Tmin=(100,'K'), Tmax=(940.446,'K')), NASAPolynomial(coeffs=[5.20299,0.00782442,-2.12683e-06,3.79691e-10,-2.94671e-14,3936.28,-6.62411], Tmin=(940.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]CCC=O(631)',
    structure = SMILES('[CH2]CCC=O'),
    E0 = (-21.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.00471129,'amu*angstrom^2'), symmetry=1, barrier=(4.90546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213589,'amu*angstrom^2'), symmetry=1, barrier=(4.91083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213759,'amu*angstrom^2'), symmetry=1, barrier=(4.91474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49762,0.0360988,-1.83384e-05,4.20693e-09,-3.76081e-13,-2578.43,18.7989], Tmin=(100,'K'), Tmax=(2478.97,'K')), NASAPolynomial(coeffs=[13.6768,0.0180607,-7.42392e-06,1.27176e-09,-8.00786e-14,-8121.12,-45.2836], Tmin=(2478.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ)"""),
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
    label = '[CH2]CC=CO(5592)',
    structure = SMILES('[CH2]CC=CO'),
    E0 = (-20.0681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.793723,'amu*angstrom^2'), symmetry=1, barrier=(18.2492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793058,'amu*angstrom^2'), symmetry=1, barrier=(18.234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792912,'amu*angstrom^2'), symmetry=1, barrier=(18.2306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64875,0.0367893,1.58547e-05,-5.77079e-08,2.81478e-11,-2315.12,20.1635], Tmin=(100,'K'), Tmax=(919.905,'K')), NASAPolynomial(coeffs=[17.227,0.00789685,-3.75915e-07,-4.01283e-11,7.62603e-17,-6824.86,-62.6259], Tmin=(919.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ)"""),
)

species(
    label = 'C=CC[CH]O(5522)',
    structure = SMILES('C=CC[CH]O'),
    E0 = (3.94862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,402.093,402.103],'cm^-1')),
        HinderedRotor(inertia=(0.0843957,'amu*angstrom^2'), symmetry=1, barrier=(9.68052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0843739,'amu*angstrom^2'), symmetry=1, barrier=(9.68377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0843046,'amu*angstrom^2'), symmetry=1, barrier=(9.68329,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8155,0.0468697,-3.90514e-05,1.80064e-08,-3.46198e-12,554.491,20.1715], Tmin=(100,'K'), Tmax=(1220.63,'K')), NASAPolynomial(coeffs=[9.24963,0.0225077,-9.11315e-06,1.65488e-09,-1.12931e-13,-1260.35,-17.1759], Tmin=(1220.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.94862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH2](22)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (313.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1607.46,1610.79,1611.1,1611.52,1612.1],'cm^-1')),
        HinderedRotor(inertia=(0.00575481,'amu*angstrom^2'), symmetry=1, barrier=(10.5915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86542,-0.00370109,4.71841e-05,-6.14996e-08,2.54414e-11,37752.3,8.63166], Tmin=(100,'K'), Tmax=(845.142,'K')), NASAPolynomial(coeffs=[5.89967,0.00441356,1.29138e-06,-4.58004e-10,3.67881e-14,36774.8,-4.58889], Tmin=(845.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ)"""),
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
    label = '[CH2]CCC[O](474)',
    structure = SMILES('[CH2]CCC[O]'),
    E0 = (126.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,184.707,1625.18,1627.48],'cm^-1')),
        HinderedRotor(inertia=(0.141284,'amu*angstrom^2'), symmetry=1, barrier=(3.98079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174723,'amu*angstrom^2'), symmetry=1, barrier=(4.0306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115059,'amu*angstrom^2'), symmetry=1, barrier=(4.03507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84954,0.0364831,-1.59529e-05,2.65371e-09,-1.22872e-13,15283.1,19.0407], Tmin=(100,'K'), Tmax=(2487.27,'K')), NASAPolynomial(coeffs=[21.9254,0.0123781,-5.3795e-06,8.82102e-10,-5.1587e-14,3760.73,-94.4578], Tmin=(2487.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[CH]CO(5593)',
    structure = SMILES('[CH2]C[CH]CO'),
    E0 = (101.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,1091.17,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00412433,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00387933,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07968,'amu*angstrom^2'), symmetry=1, barrier=(31.9326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190022,'amu*angstrom^2'), symmetry=1, barrier=(5.75876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74405,0.0444329,-2.91406e-05,1.00259e-08,-1.42473e-12,12234.9,24.6956], Tmin=(100,'K'), Tmax=(1611.09,'K')), NASAPolynomial(coeffs=[10.6959,0.0222076,-8.44795e-06,1.46338e-09,-9.60524e-14,9350.48,-22.761], Tmin=(1611.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = 'C[CH]C[CH]O(473)',
    structure = SMILES('C[CH]C[CH]O'),
    E0 = (70.6126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,188.807,2010.87],'cm^-1')),
        HinderedRotor(inertia=(0.19161,'amu*angstrom^2'), symmetry=1, barrier=(4.84139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00471741,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192072,'amu*angstrom^2'), symmetry=1, barrier=(4.83463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193906,'amu*angstrom^2'), symmetry=1, barrier=(4.83819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70787,0.0548555,-6.82338e-05,5.69033e-08,-1.97223e-11,8570.99,22.82], Tmin=(100,'K'), Tmax=(839.454,'K')), NASAPolynomial(coeffs=[4.72094,0.0325242,-1.4082e-05,2.58196e-09,-1.7474e-13,8346.08,10.4844], Tmin=(839.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.6126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCO(477)',
    structure = SMILES('[CH2][CH]CCO'),
    E0 = (95.5611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,180,1701.91],'cm^-1')),
        HinderedRotor(inertia=(0.172988,'amu*angstrom^2'), symmetry=1, barrier=(3.97734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00191122,'amu*angstrom^2'), symmetry=1, barrier=(3.94009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172301,'amu*angstrom^2'), symmetry=1, barrier=(3.96153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171809,'amu*angstrom^2'), symmetry=1, barrier=(3.95022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02007,0.046767,-4.71413e-05,3.64194e-08,-1.28128e-11,11561.6,23.9382], Tmin=(100,'K'), Tmax=(791.071,'K')), NASAPolynomial(coeffs=[3.88082,0.0334588,-1.45129e-05,2.69099e-09,-1.84481e-13,11389.2,16.1684], Tmin=(791.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.5611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH][CH]O(5594)',
    structure = SMILES('CC[CH][CH]O'),
    E0 = (76.0684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7022,0.0494489,-4.00593e-05,1.81864e-08,-3.47285e-12,9232.49,22.6002], Tmin=(100,'K'), Tmax=(1221.73,'K')), NASAPolynomial(coeffs=[9.16355,0.02502,-1.00663e-05,1.81995e-09,-1.23814e-13,7409.34,-14.8908], Tmin=(1221.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.0684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = 'CCC[CH][O](475)',
    structure = SMILES('CCC[CH][O]'),
    E0 = (101.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,248.908,249.098,1980.82],'cm^-1')),
        HinderedRotor(inertia=(0.111602,'amu*angstrom^2'), symmetry=1, barrier=(4.9132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111515,'amu*angstrom^2'), symmetry=1, barrier=(4.91129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1114,'amu*angstrom^2'), symmetry=1, barrier=(4.91059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85884,0.0520583,-6.17141e-05,5.44759e-08,-2.07038e-11,12324.6,20.4028], Tmin=(100,'K'), Tmax=(792.673,'K')), NASAPolynomial(coeffs=[3.11835,0.0370918,-1.70982e-05,3.24826e-09,-2.25146e-13,12395.4,16.3254], Tmin=(792.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]CC[CH][O](628)',
    structure = SMILES('[CH2]CC[CH][O]'),
    E0 = (307.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1530.7,1530.74],'cm^-1')),
        HinderedRotor(inertia=(0.189137,'amu*angstrom^2'), symmetry=1, barrier=(4.34863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188787,'amu*angstrom^2'), symmetry=1, barrier=(4.34058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188938,'amu*angstrom^2'), symmetry=1, barrier=(4.34406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3491.32,'J/mol'), sigma=(6.14151,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.34 K, Pc=34.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86564,0.0531628,-7.34949e-05,6.84292e-08,-2.58146e-11,37008.6,22.1639], Tmin=(100,'K'), Tmax=(823.103,'K')), NASAPolynomial(coeffs=[3.05802,0.0345392,-1.61765e-05,3.06866e-09,-2.1142e-13,37246.9,19.2834], Tmin=(823.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[CH][CH]O(746)',
    structure = SMILES('[CH2]C[CH][CH]O'),
    E0 = (281.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,180,2139.57],'cm^-1')),
        HinderedRotor(inertia=(0.010134,'amu*angstrom^2'), symmetry=1, barrier=(2.42379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105418,'amu*angstrom^2'), symmetry=1, barrier=(2.42377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0618415,'amu*angstrom^2'), symmetry=1, barrier=(14.7909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105418,'amu*angstrom^2'), symmetry=1, barrier=(2.42377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80072,0.0494004,-4.75588e-05,2.64142e-08,-6.14888e-12,33912.7,24.0379], Tmin=(100,'K'), Tmax=(1020.64,'K')), NASAPolynomial(coeffs=[8.33776,0.0237815,-9.90812e-06,1.82179e-09,-1.25207e-13,32578.3,-7.63322], Tmin=(1020.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]O(747)',
    structure = SMILES('[CH2][CH]C[CH]O'),
    E0 = (275.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,1553.81,1556.01],'cm^-1')),
        HinderedRotor(inertia=(0.181251,'amu*angstrom^2'), symmetry=1, barrier=(4.16731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181288,'amu*angstrom^2'), symmetry=1, barrier=(4.16816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180401,'amu*angstrom^2'), symmetry=1, barrier=(4.14776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181288,'amu*angstrom^2'), symmetry=1, barrier=(4.16817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72539,0.0558307,-7.95632e-05,7.03156e-08,-2.46489e-11,33254.5,24.5429], Tmin=(100,'K'), Tmax=(863.575,'K')), NASAPolynomial(coeffs=[4.59712,0.0300846,-1.32276e-05,2.41864e-09,-1.62389e-13,33222.5,13.7964], Tmin=(863.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = 'CCC=CO(5595)',
    structure = SMILES('CCC=CO'),
    E0 = (-225.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60692,0.036112,2.60996e-05,-6.96777e-08,3.24567e-11,-26997.6,18.527], Tmin=(100,'K'), Tmax=(923.192,'K')), NASAPolynomial(coeffs=[17.4571,0.0101509,-1.12125e-06,9.70601e-11,-1.0085e-14,-31744.4,-66.5329], Tmin=(923.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'C=CCCO(483)',
    structure = SMILES('C=CCCO'),
    E0 = (-176.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86161,0.0407202,-1.67183e-05,-3.01938e-09,3.01657e-12,-21127.6,20.4595], Tmin=(100,'K'), Tmax=(1102.33,'K')), NASAPolynomial(coeffs=[9.38416,0.0244476,-9.57665e-06,1.73403e-09,-1.19058e-13,-23455.9,-19.6034], Tmin=(1102.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CCCC=O(482)',
    structure = SMILES('CCCC=O'),
    E0 = (-227.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11642,0.0389099,-1.83641e-05,3.46285e-09,-1.58355e-13,-27244.3,18.4213], Tmin=(100,'K'), Tmax=(1843.53,'K')), NASAPolynomial(coeffs=[12.0049,0.0224534,-9.04166e-06,1.56257e-09,-1.00129e-13,-31739.7,-37.6376], Tmin=(1843.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
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
    label = 'OC1CCC1(5498)',
    structure = SMILES('OC1CCC1'),
    E0 = (-158.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.534,0.0153602,6.76093e-05,-9.65706e-08,3.73972e-11,-19006.8,16.7916], Tmin=(100,'K'), Tmax=(960.452,'K')), NASAPolynomial(coeffs=[11.5578,0.0203926,-6.80295e-06,1.27585e-09,-9.60547e-14,-22705.7,-36.611], Tmin=(960.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane)"""),
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
    label = '[CH]CC[CH2](104)',
    structure = SMILES('[CH]CC[CH2]'),
    E0 = (509.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,299.704,1216.25,1216.25,1216.25],'cm^-1')),
        HinderedRotor(inertia=(0.00187694,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.080812,'amu*angstrom^2'), symmetry=1, barrier=(5.15099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0808061,'amu*angstrom^2'), symmetry=1, barrier=(5.151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28697,0.0322907,-1.00182e-05,-5.96304e-09,3.57498e-12,61303.7,18.5453], Tmin=(100,'K'), Tmax=(1088.51,'K')), NASAPolynomial(coeffs=[8.15765,0.0208543,-8.22745e-06,1.49561e-09,-1.02999e-13,59425.1,-13.0337], Tmin=(1088.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C[CH2](21)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (289.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.108554,'amu*angstrom^2'), symmetry=1, barrier=(2.49587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00203216,'amu*angstrom^2'), symmetry=1, barrier=(2.52474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1302,0.0139357,1.71989e-05,-2.69382e-08,9.94977e-12,34911,13.3619], Tmin=(100,'K'), Tmax=(1011.02,'K')), NASAPolynomial(coeffs=[5.48778,0.0173102,-6.65298e-06,1.21647e-09,-8.50352e-14,33785.2,-1.24835], Tmin=(1011.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2]CC[C]O(5596)',
    structure = SMILES('[CH2]CC[C]O'),
    E0 = (362.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3100,440,815,1455,1000,180,180,1000.69],'cm^-1')),
        HinderedRotor(inertia=(0.173824,'amu*angstrom^2'), symmetry=1, barrier=(3.99655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00153523,'amu*angstrom^2'), symmetry=1, barrier=(17.4311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758223,'amu*angstrom^2'), symmetry=1, barrier=(17.433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758179,'amu*angstrom^2'), symmetry=1, barrier=(17.432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67734,0.0516592,-4.95839e-05,2.61134e-08,-5.67521e-12,43634.1,21.166], Tmin=(100,'K'), Tmax=(1094.73,'K')), NASAPolynomial(coeffs=[9.6475,0.0225377,-9.68228e-06,1.8146e-09,-1.26261e-13,41889.1,-18.0069], Tmin=(1094.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CH2_triplet) + radical(RCCJ)"""),
)

species(
    label = '[CH]CC[CH]O(5597)',
    structure = SMILES('[CH]CC[CH]O'),
    E0 = (324.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52151,0.0568352,-6.47113e-05,4.24086e-08,-1.14623e-11,39101.3,21.764], Tmin=(100,'K'), Tmax=(892.96,'K')), NASAPolynomial(coeffs=[8.72495,0.024568,-1.05094e-05,1.94305e-09,-1.33455e-13,37814.8,-12.173], Tmin=(892.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCsJOH)"""),
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
    E0 = (81.4127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (225.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (198.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (225.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (178.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (198.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (215.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (209.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (233.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (212.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (199.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (161.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (518.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (449.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (493.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (487.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (144.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (144.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (106.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (258.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (89.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (542.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (530.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (520.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (573.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (536.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['C=C(18)', 'C=CO(576)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]CCC=O(631)', 'H(8)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]CC=CO(5592)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=CC[CH]O(5522)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CO(576)', '[CH2][CH2](22)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0133609,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(18)', '[CH2][CH]O(578)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00310793,'m^3/(mol*s)'), n=2.49201, Ea=(21.4614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CCC[O](474)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C[CH]CO(5593)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['C[CH]C[CH]O(473)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]CCO(477)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['CC[CH][CH]O(5594)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC[CH][O](475)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC[CH][O](628)', 'H(8)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]O(578)', '[CH2][CH2](22)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C[CH][CH]O(746)', 'H(8)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C[CH]O(747)', 'H(8)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['CCC=CO(5595)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['C=CCCO(483)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['CCCC=O(482)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([CH2])O(3672)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CC[CH]O(5497)'],
    products = ['OC1CCC1(5498)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OH(D)(132)', '[CH]CC[CH2](104)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C[CH2](21)', '[CH]O(5471)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C[CH]O(426)', 'CH2(T)(28)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH2]CC[C]O(5596)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[CH]CC[CH]O(5597)'],
    products = ['[CH2]CC[CH]O(5497)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '1554',
    isomers = [
        '[CH2]CC[CH]O(5497)',
    ],
    reactants = [
        ('C=C(18)', 'C=CO(576)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '1554',
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


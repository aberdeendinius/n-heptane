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
    label = 'C=C1C(O)C1[C]=C[O](29165)',
    structure = SMILES('C=C1C(O)C1[C]=C[O]'),
    E0 = (256.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359748,0.0674071,-4.0847e-05,-6.63902e-09,1.05067e-11,30974.5,26.5042], Tmin=(100,'K'), Tmax=(955.419,'K')), NASAPolynomial(coeffs=[20.3419,0.0138727,-4.09241e-06,7.1516e-10,-5.27723e-14,25781.3,-76.182], Tmin=(955.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = 'C=C=C(O)C=C=C[O](29166)',
    structure = SMILES('C=C=C(O)C=C=C[O]'),
    E0 = (93.8557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33917,'amu*angstrom^2'), symmetry=1, barrier=(30.7901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32156,'amu*angstrom^2'), symmetry=1, barrier=(30.3854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20394,0.0906032,-0.000101595,5.21962e-08,-9.83046e-12,11496.1,28.9126], Tmin=(100,'K'), Tmax=(1518.86,'K')), NASAPolynomial(coeffs=[26.2571,0.00174837,2.48624e-06,-6.54997e-10,4.84411e-14,5061.38,-108.77], Tmin=(1518.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.8557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C#CC(O)C=C=C[O](29167)',
    structure = SMILES('C#CC(O)C=C=C[O]'),
    E0 = (150.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,2175,525,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2813,'amu*angstrom^2'), symmetry=1, barrier=(29.4597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27919,'amu*angstrom^2'), symmetry=1, barrier=(29.411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28249,'amu*angstrom^2'), symmetry=1, barrier=(29.4869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.294959,0.0786559,-8.34689e-05,4.25871e-08,-8.17087e-12,18222.9,31.0664], Tmin=(100,'K'), Tmax=(1431.61,'K')), NASAPolynomial(coeffs=[21.4792,0.00905502,-1.3616e-06,7.60992e-11,-5.9786e-16,12886.4,-78.6575], Tmin=(1431.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)C=C=C=O(29168)',
    structure = SMILES('C=[C]C(O)C=C=C=O'),
    E0 = (266.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.222272,'amu*angstrom^2'), symmetry=1, barrier=(5.11047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000451111,'amu*angstrom^2'), symmetry=1, barrier=(5.12194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222687,'amu*angstrom^2'), symmetry=1, barrier=(5.12001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7119,0.048172,-4.07473e-05,1.75897e-08,-3.06806e-12,32088.9,13.9195], Tmin=(100,'K'), Tmax=(1354.55,'K')), NASAPolynomial(coeffs=[11.6718,0.0187602,-8.17715e-06,1.55964e-09,-1.09496e-13,29390.6,-37.1538], Tmin=(1354.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
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
    label = 'C=C=CC=C=C[O](22623)',
    structure = SMILES('C=C=CC=C=C[O]'),
    E0 = (307.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44486,'amu*angstrom^2'), symmetry=1, barrier=(33.2202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.822392,0.0563928,-2.37759e-05,-2.27834e-08,1.65251e-11,37140.4,22.3892], Tmin=(100,'K'), Tmax=(935.614,'K')), NASAPolynomial(coeffs=[20.0281,0.00835898,-1.398e-06,1.98487e-10,-1.72302e-14,32055.2,-76.9595], Tmin=(935.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
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
    label = 'C=[C]C(O)C=C=[C]O(29169)',
    structure = SMILES('C=[C]C(O)C=C=[C]O'),
    E0 = (320.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.709985,'amu*angstrom^2'), symmetry=1, barrier=(16.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.70844,'amu*angstrom^2'), symmetry=1, barrier=(16.2884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710849,'amu*angstrom^2'), symmetry=1, barrier=(16.3438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708879,'amu*angstrom^2'), symmetry=1, barrier=(16.2985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145577,0.0814888,-9.31367e-05,5.34866e-08,-1.19848e-11,38663.5,34.0579], Tmin=(100,'K'), Tmax=(1096.37,'K')), NASAPolynomial(coeffs=[17.3312,0.0187892,-7.35463e-06,1.32567e-09,-9.08184e-14,34895.1,-50.4339], Tmin=(1096.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CJO)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4502.04,'J/mol'), sigma=(7.09037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=703.21 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438181,0.0654159,-3.65856e-05,-7.54109e-09,9.51967e-12,25941.3,30.6225], Tmin=(100,'K'), Tmax=(987.953,'K')), NASAPolynomial(coeffs=[19.5037,0.016784,-6.11062e-06,1.15527e-09,-8.53236e-14,20780.3,-68.1804], Tmin=(987.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)OJ)"""),
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
    label = 'C=[C]C(O)C#C[CH]O(29170)',
    structure = SMILES('C=[C]C(O)C#C[CH]O'),
    E0 = (299.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,2100,2250,500,550,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595624,0.079287,-9.59035e-05,5.43683e-08,-7.93505e-12,36163,32.9689], Tmin=(100,'K'), Tmax=(686.228,'K')), NASAPolynomial(coeffs=[12.4453,0.0249518,-9.34531e-06,1.57022e-09,-1.00467e-13,34189.7,-22.265], Tmin=(686.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S) + radical(CCsJOH)"""),
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
    label = 'C=[C]C(O)=C[C]=CO(29171)',
    structure = SMILES('C=[C]C(O)=C[C]=CO'),
    E0 = (120.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01139,0.106903,-0.000131028,7.30527e-08,-1.47142e-11,14791.1,30.4881], Tmin=(100,'K'), Tmax=(1455.22,'K')), NASAPolynomial(coeffs=[28.323,-0.00133086,6.15323e-06,-1.52859e-09,1.14632e-13,8593.89,-118.197], Tmin=(1455.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = '[CH]=[C]C(O)C=C=CO(29172)',
    structure = SMILES('[CH]=[C]C(O)C=C=CO'),
    E0 = (327.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,540,610,2055,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.900574,'amu*angstrom^2'), symmetry=1, barrier=(20.706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900856,'amu*angstrom^2'), symmetry=1, barrier=(20.7125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900782,'amu*angstrom^2'), symmetry=1, barrier=(20.7108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900791,'amu*angstrom^2'), symmetry=1, barrier=(20.711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.55312,0.0881521,-0.000100057,5.43463e-08,-1.11972e-11,39580.2,33.8391], Tmin=(100,'K'), Tmax=(1263.83,'K')), NASAPolynomial(coeffs=[22.9996,0.0096742,-2.24477e-06,2.87481e-10,-1.65282e-14,33941.1,-84.0607], Tmin=(1263.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C=C[C]=C[O](22614)',
    structure = SMILES('C=[C]C=C[C]=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3934.94,'J/mol'), sigma=(6.27619,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.63 K, Pc=36.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = 'C=[C]C(O)=C[C]=C[O](29173)',
    structure = SMILES('C=[C]C(O)=C[C]=C[O]'),
    E0 = (262.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.46364,'amu*angstrom^2'), symmetry=1, barrier=(33.652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4765,'amu*angstrom^2'), symmetry=1, barrier=(33.9476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47159,'amu*angstrom^2'), symmetry=1, barrier=(33.8346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35286,0.0969676,-0.000116984,6.50287e-08,-1.31748e-11,31777.4,29.5381], Tmin=(100,'K'), Tmax=(1427.58,'K')), NASAPolynomial(coeffs=[25.7624,0.00141754,3.98105e-06,-1.06588e-09,8.18072e-14,26030.2,-103.944], Tmin=(1427.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C(O)[C]=C=C[O](29174)',
    structure = SMILES('C=[C]C(O)[C]=C=C[O]'),
    E0 = (459.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.837507,'amu*angstrom^2'), symmetry=1, barrier=(19.2559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836238,'amu*angstrom^2'), symmetry=1, barrier=(19.2267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83727,'amu*angstrom^2'), symmetry=1, barrier=(19.2505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359182,0.0768218,-8.57087e-05,4.76085e-08,-1.03349e-11,55441,32.0469], Tmin=(100,'K'), Tmax=(1128.58,'K')), NASAPolynomial(coeffs=[16.9676,0.0179568,-7.47108e-06,1.39256e-09,-9.72578e-14,51692.2,-50.0881], Tmin=(1128.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C(O)C=C=C[O](29175)',
    structure = SMILES('[CH]=[C]C(O)C=C=C[O]'),
    E0 = (469.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.9137,'amu*angstrom^2'), symmetry=1, barrier=(21.0078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910768,'amu*angstrom^2'), symmetry=1, barrier=(20.9404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913057,'amu*angstrom^2'), symmetry=1, barrier=(20.993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136001,0.077892,-8.50491e-05,4.53206e-08,-9.33104e-12,56565.3,32.7771], Tmin=(100,'K'), Tmax=(1196.58,'K')), NASAPolynomial(coeffs=[19.1263,0.0144103,-5.47025e-06,9.83854e-10,-6.78299e-14,52020.6,-62.2481], Tmin=(1196.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(O)[CH][C]=C=O(29176)',
    structure = SMILES('C=[C]C(O)[CH][C]=C=O'),
    E0 = (315.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2120,512.5,787.5,1670,1700,300,440,1380,1390,370,380,2900,435,366.607,366.971,367.09],'cm^-1')),
        HinderedRotor(inertia=(0.00125279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347104,'amu*angstrom^2'), symmetry=1, barrier=(33.2524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348507,'amu*angstrom^2'), symmetry=1, barrier=(33.2618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348371,'amu*angstrom^2'), symmetry=1, barrier=(33.2544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658241,0.0626682,-4.19912e-05,3.78531e-09,4.27433e-12,38135.5,32.0207], Tmin=(100,'K'), Tmax=(1023.71,'K')), NASAPolynomial(coeffs=[17.6911,0.0174429,-6.97517e-06,1.33332e-09,-9.68514e-14,33530.6,-56.0103], Tmin=(1023.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C=[C]C(O)C=C1[CH]O1(29177)',
    structure = SMILES('C=[C]C(O)C=C1[CH]O1'),
    E0 = (248.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661317,0.0528515,9.77035e-06,-6.53454e-08,3.26753e-11,29989.3,29.0679], Tmin=(100,'K'), Tmax=(943.906,'K')), NASAPolynomial(coeffs=[23.4163,0.00802748,-1.00588e-06,1.86734e-10,-2.21261e-14,23394.7,-91.5759], Tmin=(943.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = 'C=C1C([CH]C1O)=C[O](29117)',
    structure = SMILES('C=C1C([CH]C1O)=C[O]'),
    E0 = (16.1651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991248,0.0356235,7.61164e-05,-1.42599e-07,6.2003e-11,2080.7,26.1035], Tmin=(100,'K'), Tmax=(934.16,'K')), NASAPolynomial(coeffs=[26.7027,0.00358546,2.2239e-06,-4.17934e-10,1.47642e-14,-6128.82,-114.417], Tmin=(934.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.1651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=CCJC(O)C=C) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)C1[C]=CO1(29178)',
    structure = SMILES('C=[C]C(O)C1[C]=CO1'),
    E0 = (345.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444859,0.055877,9.13897e-06,-7.04716e-08,3.5949e-11,41714.8,29.7146], Tmin=(100,'K'), Tmax=(937.021,'K')), NASAPolynomial(coeffs=[25.7534,0.00442104,9.32431e-07,-1.88778e-10,3.45678e-15,34487.9,-103.993], Tmin=(937.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C([O])[C]=CC1O(29179)',
    structure = SMILES('C=C1C([O])[C]=CC1O'),
    E0 = (247.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31182,0.0472782,-5.96926e-06,-2.35604e-08,1.11928e-11,29882.9,28.0135], Tmin=(100,'K'), Tmax=(1063.32,'K')), NASAPolynomial(coeffs=[13.9202,0.0253497,-1.10103e-05,2.15548e-09,-1.56437e-13,25759.9,-40.3677], Tmin=(1063.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C=C=C(O)C=C=CO(29180)',
    structure = SMILES('C=C=C(O)C=C=CO'),
    E0 = (-47.6069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87083,0.100625,-0.000115892,6.04879e-08,-1.14611e-11,-5489.9,29.8934], Tmin=(100,'K'), Tmax=(1534.46,'K')), NASAPolynomial(coeffs=[28.7005,-0.000842766,4.58285e-06,-1.10222e-09,8.01164e-14,-12308.5,-122.333], Tmin=(1534.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.6069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH2][C]=C(O)C[C]=C[O](29181)',
    structure = SMILES('[CH2][C]=C(O)C[C]=C[O]'),
    E0 = (398.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,229.855,229.857,229.862],'cm^-1')),
        HinderedRotor(inertia=(0.53915,'amu*angstrom^2'), symmetry=1, barrier=(20.2173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539152,'amu*angstrom^2'), symmetry=1, barrier=(20.2172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539166,'amu*angstrom^2'), symmetry=1, barrier=(20.2173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539219,'amu*angstrom^2'), symmetry=1, barrier=(20.2173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.888063,0.0869799,-9.26608e-05,4.65159e-08,-8.69591e-12,48064.1,34.738], Tmin=(100,'K'), Tmax=(1494.61,'K')), NASAPolynomial(coeffs=[24.5084,0.00692459,-1.86113e-07,-1.43112e-10,1.37008e-14,41822.6,-93.4749], Tmin=(1494.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(O)[C]C=C[O](29182)',
    structure = SMILES('C=[C]C(O)[C]C=C[O]'),
    E0 = (509.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.1053,0.0749129,-6.22395e-05,1.57627e-08,2.35022e-12,61443,33.9352], Tmin=(100,'K'), Tmax=(976.33,'K')), NASAPolynomial(coeffs=[20.5496,0.0143382,-4.7953e-06,8.60996e-10,-6.2172e-14,56345.9,-69.866], Tmin=(976.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]C(O)[CH]C=[C][O](29183)',
    structure = SMILES('C=[C]C(O)[CH]C=[C][O]'),
    E0 = (371.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1670,1700,300,440,3010,987.5,1337.5,450,1655,561.504,561.505,561.505,561.505],'cm^-1')),
        HinderedRotor(inertia=(0.0965178,'amu*angstrom^2'), symmetry=1, barrier=(21.5944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965182,'amu*angstrom^2'), symmetry=1, barrier=(21.5944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965175,'amu*angstrom^2'), symmetry=1, barrier=(21.5944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965178,'amu*angstrom^2'), symmetry=1, barrier=(21.5944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631031,0.0593351,-1.92277e-05,-2.60322e-08,1.62529e-11,44790.8,36.0704], Tmin=(100,'K'), Tmax=(973.107,'K')), NASAPolynomial(coeffs=[19.6205,0.0158409,-5.46037e-06,1.03555e-09,-7.81809e-14,39458.6,-63.4334], Tmin=(973.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
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
    label = '[CH2][C]=C(O)[CH]C=C[O](29184)',
    structure = SMILES('[CH2][C]=C(O)[CH]C=C[O]'),
    E0 = (277.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3434,'amu*angstrom^2'), symmetry=1, barrier=(30.8874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34294,'amu*angstrom^2'), symmetry=1, barrier=(30.8768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3436,'amu*angstrom^2'), symmetry=1, barrier=(30.892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34138,'amu*angstrom^2'), symmetry=1, barrier=(30.8411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287268,0.0772045,-4.42774e-05,-1.88434e-08,1.86516e-11,33496.8,29.2984], Tmin=(100,'K'), Tmax=(929.823,'K')), NASAPolynomial(coeffs=[26.4798,0.00583349,2.3606e-07,-1.23425e-10,4.12498e-15,26626.6,-108.066], Tmin=(929.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJCO) + radical(Cds_S) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = '[CH2][C]=C(O)C=[C]C[O](29185)',
    structure = SMILES('[CH2][C]=C(O)C=[C]C[O]'),
    E0 = (474.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,262.538,262.538,3077.78],'cm^-1')),
        HinderedRotor(inertia=(0.280228,'amu*angstrom^2'), symmetry=1, barrier=(13.7064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48047,'amu*angstrom^2'), symmetry=1, barrier=(72.4126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28023,'amu*angstrom^2'), symmetry=1, barrier=(13.7064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48048,'amu*angstrom^2'), symmetry=1, barrier=(72.4126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434905,0.0821215,-0.000108158,7.9658e-08,-2.36667e-11,57168.3,30.7726], Tmin=(100,'K'), Tmax=(822.272,'K')), NASAPolynomial(coeffs=[11.1841,0.0298338,-1.27783e-05,2.33125e-09,-1.57668e-13,55400.5,-18.9831], Tmin=(822.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ) + radical(Cds_S) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[C]=C(584)',
    structure = SMILES('[C]=C'),
    E0 = (600.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94093,-0.00117598,1.80376e-05,-2.01208e-08,6.96659e-12,72197.9,5.25681], Tmin=(100,'K'), Tmax=(976.125,'K')), NASAPolynomial(coeffs=[3.93016,0.00536132,-1.98619e-06,3.69549e-10,-2.66221e-14,71890.7,3.724], Tmin=(976.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[O]C=[C]C=CO(23185)',
    structure = SMILES('[O]C=[C]C=CO'),
    E0 = (16.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44835,'amu*angstrom^2'), symmetry=1, barrier=(33.3004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45066,'amu*angstrom^2'), symmetry=1, barrier=(33.3536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93068,0.0480724,2.96789e-06,-6.71246e-08,3.83779e-11,2148.87,19.742], Tmin=(100,'K'), Tmax=(884.221,'K')), NASAPolynomial(coeffs=[26.0788,-0.010819,9.78516e-06,-2.08098e-09,1.44498e-13,-4443.54,-110.619], Tmin=(884.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = '[CH]=C=CC(O)[C]=C(27788)',
    structure = SMILES('[CH]=C=CC(O)[C]=C'),
    E0 = (443.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.694272,'amu*angstrom^2'), symmetry=1, barrier=(15.9627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694783,'amu*angstrom^2'), symmetry=1, barrier=(15.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694709,'amu*angstrom^2'), symmetry=1, barrier=(15.9727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06512,0.0640079,-6.62022e-05,3.69969e-08,-8.35491e-12,53484.7,28.3977], Tmin=(100,'K'), Tmax=(1069.83,'K')), NASAPolynomial(coeffs=[11.8516,0.0236781,-9.65573e-06,1.75967e-09,-1.20544e-13,51176.7,-24.3687], Tmin=(1069.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)C#CC=O(29186)',
    structure = SMILES('C=[C]C(O)C#CC=O'),
    E0 = (184.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.063787,'amu*angstrom^2'), symmetry=1, barrier=(15.1866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660218,'amu*angstrom^2'), symmetry=1, barrier=(15.1797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06767,'amu*angstrom^2'), symmetry=1, barrier=(87.5728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0638019,'amu*angstrom^2'), symmetry=1, barrier=(15.1827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05077,0.067131,-7.26999e-05,4.30602e-08,-1.04439e-11,22324.3,29.4132], Tmin=(100,'K'), Tmax=(990.932,'K')), NASAPolynomial(coeffs=[10.9663,0.0271067,-1.2115e-05,2.30136e-09,-1.611e-13,20359.2,-18.3333], Tmin=(990.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[C]=CC=O(29187)',
    structure = SMILES('C=[C]C(O)[C]=CC=O'),
    E0 = (263.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1670,1700,300,440,3010,987.5,1337.5,450,1655,288.963,288.968],'cm^-1')),
        HinderedRotor(inertia=(0.184783,'amu*angstrom^2'), symmetry=1, barrier=(10.9508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184796,'amu*angstrom^2'), symmetry=1, barrier=(10.9508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184806,'amu*angstrom^2'), symmetry=1, barrier=(10.9508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184787,'amu*angstrom^2'), symmetry=1, barrier=(10.9508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827141,0.0769833,-0.000102592,8.25109e-08,-2.82446e-11,31752.6,31.3266], Tmin=(100,'K'), Tmax=(701.858,'K')), NASAPolynomial(coeffs=[7.61475,0.038305,-1.99408e-05,4.01515e-09,-2.88536e-13,30799.7,0.982332], Tmin=(701.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(O)[CH]C=C=O(29188)',
    structure = SMILES('C=[C]C(O)[CH]C=C=O'),
    E0 = (113.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676242,0.0615778,-3.45163e-05,-3.45961e-09,6.51467e-12,13806.7,32.2474], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[17.2717,0.0201833,-8.09612e-06,1.54624e-09,-1.12e-13,9175.63,-54.2176], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=[C]C(O)=CC=C[O](29189)',
    structure = SMILES('C=[C]C(O)=CC=C[O]'),
    E0 = (63.4657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78284,0.0960774,-0.000106196,5.39809e-08,-9.95303e-12,7868.17,31.5079], Tmin=(100,'K'), Tmax=(1594.73,'K')), NASAPolynomial(coeffs=[26.615,0.00199293,3.79709e-06,-9.87976e-10,7.30897e-14,1717.02,-109.637], Tmin=(1594.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.4657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = 'C#CC(O)[CH]C=C[O](25783)',
    structure = SMILES('C#CC(O)[CH]C=C[O]'),
    E0 = (106.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0238042,0.0737386,-4.33502e-05,-1.4298e-08,1.59966e-11,12959.1,28.8364], Tmin=(100,'K'), Tmax=(924.681,'K')), NASAPolynomial(coeffs=[23.6945,0.00944229,-1.18738e-06,1.03077e-10,-8.96117e-15,6935.14,-92.5885], Tmin=(924.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1OC=[C][CH]C1O(29190)',
    structure = SMILES('C=C1OC=[C][CH]C1O'),
    E0 = (79.4421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63472,0.0213343,0.000105949,-1.58239e-07,6.24929e-11,9667.56,24.0513], Tmin=(100,'K'), Tmax=(966.94,'K')), NASAPolynomial(coeffs=[21.6424,0.0158519,-5.43694e-06,1.21742e-09,-1.0566e-13,2185.35,-90.4842], Tmin=(966.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(O)C=CC=O(29191)',
    structure = SMILES('C=C=C(O)C=CC=O'),
    E0 = (-102.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0935441,0.0811927,-8.55872e-05,4.43348e-08,-8.98635e-12,-12229.1,25.1987], Tmin=(100,'K'), Tmax=(1203.51,'K')), NASAPolynomial(coeffs=[18.704,0.0193389,-8.4957e-06,1.63125e-09,-1.15765e-13,-16708.7,-68.0336], Tmin=(1203.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C1C(C=O)=CC1O(29137)',
    structure = SMILES('C=C1C(C=O)=CC1O'),
    E0 = (-72.0625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.664756,0.0556382,-5.08591e-06,-4.20168e-08,2.17695e-11,-8531.2,23.4139], Tmin=(100,'K'), Tmax=(984.516,'K')), NASAPolynomial(coeffs=[21.1206,0.0141631,-5.32978e-06,1.10346e-09,-8.7744e-14,-14576.8,-85.2019], Tmin=(984.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.0625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
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
    label = '[C]=CC(O)[C]=C(28219)',
    structure = SMILES('[C]=CC(O)[C]=C'),
    E0 = (706.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,242.434,242.437],'cm^-1')),
        HinderedRotor(inertia=(0.00286792,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284192,'amu*angstrom^2'), symmetry=1, barrier=(11.8547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284217,'amu*angstrom^2'), symmetry=1, barrier=(11.8548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80339,0.051898,-5.73104e-05,3.59044e-08,-9.43315e-12,85091.7,24.6497], Tmin=(100,'K'), Tmax=(907.323,'K')), NASAPolynomial(coeffs=[8.08798,0.0241917,-1.15059e-05,2.24882e-09,-1.59781e-13,83951.3,-5.05859], Tmin=(907.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
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
    E0 = (222.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (305.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (318.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (375.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (480.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (265.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (340.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (404.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (336.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (497.975,'kJ/mol'),
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
    E0 = (380.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (614.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (314.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (313.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (383.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (479.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (504.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (671.267,'kJ/mol'),
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
    E0 = (474.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (671.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (680.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (527.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (409.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (348.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (347.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (278.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (239.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (261.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (423.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (531.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (394.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (607.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (340.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (386.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (482.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (678.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (316.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (651.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (850.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (406.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (287.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (459.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (384.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (398.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (675.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (464.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (280.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (315.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (230.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (774.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C=CO(12571)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C1C(O)C1[C]=C[O](29165)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=C=C(O)C=C=C[O](29166)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC(O)C=C=C[O](29167)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C=[C]C(O)C=C=C=O(29168)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=C[O](8556)', 'C=C=CO(12571)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00401797,'m^3/(mol*s)'), n=2.41733, Ea=(22.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;CJ] for rate rule [Cds_Ca;CdsJ=Cdd]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(D)(132)', 'C=C=CC=C=C[O](22623)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.971254,'m^3/(mol*s)'), n=2.031, Ea=(4.66604,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Ca;OJ_pri] + [Cds-CdH_Ca;YJ] for rate rule [Cds-CdH_Ca;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['[CH2]C=C(O)C=C=C[O](25842)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC(O)C=C=C[O](25846)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(O)C=C=[C]O(29169)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
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
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.85004e+08,'s^-1'), n=1.20667, Ea=(158.922,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] + [R3H_SS_Cs;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(O)C#C[CH]O(29170)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=CC(O)[CH][C]=C=O(25848)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.08955e+06,'s^-1'), n=1.58778, Ea=(92.3497,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R5H;Cd_rad_out;XH_out] for rate rule [R5H;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=[C]C(O)=C[C]=CO(29171)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C([O])C=C=CO(25850)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;O_rad_out;XH_out] for rate rule [R6H;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(O)C=C=CO(29172)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OH(D)(132)', 'C=[C]C=C[C]=C[O](22614)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', 'C=[C]C([O])C=C=C[O](25853)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.00518e+06,'m^3/(mol*s)'), n=0.282325, Ea=(7.09479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=C[O](8556)', '[CH2][C]=CO(18753)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', 'C=[C]C(O)=C[C]=C[O](29173)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C=[C]C(O)[C]=C=C[O](29174)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', '[CH]=[C]C(O)C=C=C[O](29175)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', 'C=[C]C(O)[CH][C]=C=O(29176)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=[C]C(O)C=C1[CH]O1(29177)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C1C([CH]C1O)=C[O](29117)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=[C]C(O)C1[C]=CO1(29178)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C1C([O])[C]=CC1O(29179)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.01659e+11,'s^-1'), n=0.0800815, Ea=(56.4798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra;radadd_intra_cddouble] + [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C=C(O)C=C=CO(29180)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=CC(O)C=C=C=O(25861)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=C(O)C[C]=C[O](29181)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(O)[C]C=C[O](29182)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(O)[CH]C=[C][O](29183)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C([O])C[C]=C[O](25866)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C]=C(O)[CH]C=C[O](29184)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]C([O])[CH]C=C[O](25869)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C]=C(O)C=[C]C[O](29185)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C]C([O])C=[C]C[O](25873)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['[CH2]C(C=C=C[O])=CO(25407)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[C]=C(584)', '[O]C=[C]C=CO(23185)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O(T)(63)', '[CH]=C=CC(O)[C]=C(27788)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(8)', 'C=[C]C(O)C#CC=O(29186)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][C]=CO(18753)', 'C#CC=O(21959)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C]C(O)[C]=CC=O(29187)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=[C]C(O)[CH]C=C=O(29188)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=[C]C(O)=CC=C[O](29189)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C]C([O])C=CC=O(25883)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C#CC(O)[CH]C=C[O](25783)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.09763e+09,'s^-1'), n=1.36198, Ea=(242.913,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleDe;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C1OC=[C][CH]C1O(29190)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.19745e+11,'s^-1'), n=0.440371, Ea=(58.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C=C(O)C=CC=O(29191)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=[C]C(O)C=C=C[O](25844)'],
    products = ['C=C1C(C=O)=CC1O(29137)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=O(373)', '[C]=CC(O)[C]=C(28219)'],
    products = ['C=[C]C(O)C=C=C[O](25844)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '5260',
    isomers = [
        'C=[C]C(O)C=C=C[O](25844)',
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
    label = '5260',
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


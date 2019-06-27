species(
    label = 'C=C([O])CC=C[O](12760)',
    structure = SMILES('C=C([O])CC=C[O]'),
    E0 = (-55.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,300.227,301.17,301.318,301.517],'cm^-1')),
        HinderedRotor(inertia=(0.27198,'amu*angstrom^2'), symmetry=1, barrier=(17.397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273679,'amu*angstrom^2'), symmetry=1, barrier=(17.4024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.978133,0.0515904,-7.28955e-06,-3.97615e-08,2.25807e-11,-6532.62,26.4979], Tmin=(100,'K'), Tmax=(929.656,'K')), NASAPolynomial(coeffs=[19.5071,0.00995383,-1.56325e-06,2.01749e-10,-1.71757e-14,-11623.6,-70.3944], Tmin=(929.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C=O(598)',
    structure = SMILES('C=C=O'),
    E0 = (-59.3981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52746,0.00708371,9.17709e-06,-1.64254e-08,6.71115e-12,-7123.94,5.7438], Tmin=(100,'K'), Tmax=(956.683,'K')), NASAPolynomial(coeffs=[5.76495,0.00596559,-1.98486e-06,3.52744e-10,-2.51619e-14,-7929,-6.92178], Tmin=(956.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.3981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-(Cdd-O2d)HH)"""),
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
    label = 'C=C1CC([CH][O])O1(14462)',
    structure = SMILES('C=C1CC([CH][O])O1'),
    E0 = (154.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13516,0.0525005,-1.81305e-05,-2.29302e-08,1.57708e-11,18696.8,20.2124], Tmin=(100,'K'), Tmax=(909.979,'K')), NASAPolynomial(coeffs=[15.9988,0.0154053,-3.53564e-06,4.82556e-10,-3.12146e-14,14822.4,-56.5188], Tmin=(909.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1([O])CC=CO1(14440)',
    structure = SMILES('[CH2]C1([O])CC=CO1'),
    E0 = (65.3864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16237,0.0398192,4.19886e-05,-1.03828e-07,4.9543e-11,7987.61,18.5832], Tmin=(100,'K'), Tmax=(898.717,'K')), NASAPolynomial(coeffs=[23.4435,0.00134801,4.89203e-06,-1.16052e-09,7.9024e-14,1531.49,-100.168], Tmin=(898.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.3864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)(O)OJ) + radical(CJC(C)OC)"""),
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
    label = 'C=C([O])CC=C=O(14463)',
    structure = SMILES('C=C([O])CC=C=O'),
    E0 = (-74.4878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2120,512.5,787.5,3010,987.5,1337.5,450,1655,214.966,216.665,221.928],'cm^-1')),
        HinderedRotor(inertia=(0.236901,'amu*angstrom^2'), symmetry=1, barrier=(8.4824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241269,'amu*angstrom^2'), symmetry=1, barrier=(8.43784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18792,0.0620475,-6.91499e-05,4.13431e-08,-9.90986e-12,-8857.67,23.7413], Tmin=(100,'K'), Tmax=(1014.27,'K')), NASAPolynomial(coeffs=[11.5614,0.0211374,-8.64774e-06,1.57562e-09,-1.07816e-13,-10962,-26.4514], Tmin=(1014.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.4878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])CC=[C]O(14464)',
    structure = SMILES('C=C([O])CC=[C]O'),
    E0 = (42.9475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,297.788,297.789,297.789],'cm^-1')),
        HinderedRotor(inertia=(0.211035,'amu*angstrom^2'), symmetry=1, barrier=(13.28,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211034,'amu*angstrom^2'), symmetry=1, barrier=(13.28,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211035,'amu*angstrom^2'), symmetry=1, barrier=(13.28,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0410214,0.069785,-6.99637e-05,3.39517e-08,-6.14125e-12,5327.28,31.8857], Tmin=(100,'K'), Tmax=(1571.91,'K')), NASAPolynomial(coeffs=[19.4154,0.00876416,-7.50481e-07,-6.09036e-11,9.06892e-15,632.593,-66.2569], Tmin=(1571.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[CH]C=C[O](14465)',
    structure = SMILES('C=C(O)[CH]C=C[O]'),
    E0 = (-76.2221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400843,0.0594199,-7.08619e-06,-5.37493e-08,3.0632e-11,-9019.51,24.2686], Tmin=(100,'K'), Tmax=(924.396,'K')), NASAPolynomial(coeffs=[25.3551,0.00260988,2.06445e-06,-4.65225e-10,2.62503e-14,-15819.3,-105.985], Tmin=(924.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.2221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C(O)CC=C[O](14466)',
    structure = SMILES('[CH]=C(O)CC=C[O]'),
    E0 = (53.9573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.927047,'amu*angstrom^2'), symmetry=1, barrier=(21.3146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.930049,'amu*angstrom^2'), symmetry=1, barrier=(21.3837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928544,'amu*angstrom^2'), symmetry=1, barrier=(21.349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546559,0.058973,-1.59416e-05,-4.02783e-08,2.5137e-11,6629.69,26.9498], Tmin=(100,'K'), Tmax=(920.131,'K')), NASAPolynomial(coeffs=[23.4752,0.00378054,1.51645e-06,-3.86316e-10,2.29209e-14,527.145,-91.9919], Tmin=(920.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=C([O])C[C]=CO(14467)',
    structure = SMILES('C=C([O])C[C]=CO'),
    E0 = (41.0451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.749384,'amu*angstrom^2'), symmetry=1, barrier=(17.2298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750192,'amu*angstrom^2'), symmetry=1, barrier=(17.2484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75198,'amu*angstrom^2'), symmetry=1, barrier=(17.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.56428,0.0758759,-7.90326e-05,3.87517e-08,-6.96524e-12,5122.04,31.1132], Tmin=(100,'K'), Tmax=(1620.33,'K')), NASAPolynomial(coeffs=[21.4332,0.00536805,1.23983e-06,-4.47392e-10,3.50195e-14,120.562,-79.0647], Tmin=(1620.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.0451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C[C]=C[O](14468)',
    structure = SMILES('C=C(O)C[C]=C[O]'),
    E0 = (44.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.868046,'amu*angstrom^2'), symmetry=1, barrier=(19.9581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86621,'amu*angstrom^2'), symmetry=1, barrier=(19.9159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86788,'amu*angstrom^2'), symmetry=1, barrier=(19.9543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59056,0.0600186,-2.40038e-05,-2.83974e-08,2.0043e-11,5513.16,26.8621], Tmin=(100,'K'), Tmax=(923.197,'K')), NASAPolynomial(coeffs=[22.0147,0.00616495,1.75486e-07,-1.31694e-10,6.16194e-15,-103.359,-83.7796], Tmin=(923.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=C([O])[CH]C=CO(14469)',
    structure = SMILES('C=C([O])[CH]C=CO'),
    E0 = (-79.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347664,0.0625043,-1.82904e-05,-4.25449e-08,2.72266e-11,-9459.18,24.5509], Tmin=(100,'K'), Tmax=(911.865,'K')), NASAPolynomial(coeffs=[24.9673,0.00211818,2.72483e-06,-6.50021e-10,4.22106e-14,-15928.6,-102.807], Tmin=(911.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C[CH][C]=O(14470)',
    structure = SMILES('C=C(O)C[CH][C]=O'),
    E0 = (-10.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,1855,455,950,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,307.029,3617.44],'cm^-1')),
        HinderedRotor(inertia=(0.29235,'amu*angstrom^2'), symmetry=1, barrier=(19.5547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292297,'amu*angstrom^2'), symmetry=1, barrier=(19.5548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292299,'amu*angstrom^2'), symmetry=1, barrier=(19.5549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292346,'amu*angstrom^2'), symmetry=1, barrier=(19.5546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303109,0.0682155,-6.86573e-05,3.42762e-08,-6.51411e-12,-1078.89,28.2644], Tmin=(100,'K'), Tmax=(1426.93,'K')), NASAPolynomial(coeffs=[18.152,0.0118123,-2.67081e-06,3.19224e-10,-1.67854e-14,-5524.33,-61.92], Tmin=(1426.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C([O])CC=CO(14471)',
    structure = SMILES('[CH]=C([O])CC=CO'),
    E0 = (50.2995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.824943,'amu*angstrom^2'), symmetry=1, barrier=(18.9671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821267,'amu*angstrom^2'), symmetry=1, barrier=(18.8825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823427,'amu*angstrom^2'), symmetry=1, barrier=(18.9322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.875274,0.077837,-8.09191e-05,3.91244e-08,-6.87577e-12,6250.58,32.1695], Tmin=(100,'K'), Tmax=(1681.92,'K')), NASAPolynomial(coeffs=[22.0134,0.0040816,2.08988e-06,-6.08073e-10,4.52478e-14,1284.01,-82.0315], Tmin=(1681.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.2995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=O(601)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,672.051,672.102],'cm^-1')),
        HinderedRotor(inertia=(0.000373196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57974,0.00389613,2.17609e-05,-3.06386e-08,1.18311e-11,19367.5,10.1348], Tmin=(100,'K'), Tmax=(961.532,'K')), NASAPolynomial(coeffs=[6.4326,0.00553733,-1.87382e-06,3.59985e-10,-2.76653e-14,18194.3,-6.76404], Tmin=(961.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C[O](602)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (221.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27415,0.00479611,3.742e-05,-6.11894e-08,2.65903e-11,26726.9,9.63858], Tmin=(100,'K'), Tmax=(905.806,'K')), NASAPolynomial(coeffs=[11.9892,-0.00434473,3.96329e-06,-8.00891e-10,5.23184e-14,23944.2,-38.1893], Tmin=(905.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=C)[O](4273)',
    structure = SMILES('[CH2]C(=C)[O]'),
    E0 = (88.2866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,510.595],'cm^-1')),
        HinderedRotor(inertia=(0.0480287,'amu*angstrom^2'), symmetry=1, barrier=(8.88265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.98,'J/mol'), sigma=(5.64088,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.76 K, Pc=42.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6374,0.0235792,5.32605e-07,-2.30624e-08,1.26355e-11,10673.5,14.3058], Tmin=(100,'K'), Tmax=(894.06,'K')), NASAPolynomial(coeffs=[10.3562,0.00670937,-7.99446e-07,2.86693e-11,-3.46262e-16,8587.33,-26.0166], Tmin=(894.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.2866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])[CH]C=C[O](14472)',
    structure = SMILES('C=C([O])[CH]C=C[O]'),
    E0 = (61.5827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.4463,'amu*angstrom^2'), symmetry=1, barrier=(33.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44769,'amu*angstrom^2'), symmetry=1, barrier=(33.2852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762041,0.0555156,-1.49125e-05,-3.61221e-08,2.22965e-11,7537.78,24.4728], Tmin=(100,'K'), Tmax=(928.59,'K')), NASAPolynomial(coeffs=[21.632,0.00594283,2.28186e-08,-7.68365e-11,1.09364e-15,1923.19,-84.0284], Tmin=(928.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.5827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C([O])C[C]=C[O](14473)',
    structure = SMILES('C=C([O])C[C]=C[O]'),
    E0 = (182.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,317.063,317.388,317.474,317.483],'cm^-1')),
        HinderedRotor(inertia=(0.225894,'amu*angstrom^2'), symmetry=1, barrier=(16.1114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225335,'amu*angstrom^2'), symmetry=1, barrier=(16.1102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950636,0.0561282,-3.18828e-05,-1.06961e-08,1.16732e-11,22070.5,27.0702], Tmin=(100,'K'), Tmax=(928.985,'K')), NASAPolynomial(coeffs=[18.2946,0.0094926,-1.86303e-06,2.55947e-10,-1.89319e-14,17637.9,-61.84], Tmin=(928.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C([O])CC=C[O](14474)',
    structure = SMILES('[CH]=C([O])CC=C[O]'),
    E0 = (191.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,198.328,198.348,198.39],'cm^-1')),
        HinderedRotor(inertia=(0.645382,'amu*angstrom^2'), symmetry=1, barrier=(17.9971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644854,'amu*angstrom^2'), symmetry=1, barrier=(17.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907963,0.0550667,-2.37639e-05,-2.26507e-08,1.67981e-11,23187,27.1532], Tmin=(100,'K'), Tmax=(923.733,'K')), NASAPolynomial(coeffs=[19.7491,0.0071186,-5.2811e-07,2.76001e-12,-2.29249e-15,18271,-70.0182], Tmin=(923.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]([O])CC=C=O(10744)',
    structure = SMILES('[CH2][C]([O])CC=C=O'),
    E0 = (320.644,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3010,987.5,1337.5,450,1655,360,370,350,229.115,229.737,230.279],'cm^-1')),
        HinderedRotor(inertia=(0.00313984,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318785,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295834,'amu*angstrom^2'), symmetry=1, barrier=(11.2351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.838125,0.0744895,-0.000108964,8.81576e-08,-2.8606e-11,38673.7,26.6333], Tmin=(100,'K'), Tmax=(805.525,'K')), NASAPolynomial(coeffs=[9.77345,0.0263433,-1.22773e-05,2.31899e-09,-1.59352e-13,37356.7,-13.7817], Tmin=(805.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[O]C=CC[C]1CO1(9768)',
    structure = SMILES('[O]C=CC[C]1CO1'),
    E0 = (89.5659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,3150,900,1100,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0219926,0.0655304,-6.12443e-05,2.96917e-08,-5.34876e-12,10934.5,28.3911], Tmin=(100,'K'), Tmax=(1650.19,'K')), NASAPolynomial(coeffs=[14.6091,0.0148873,-1.3171e-06,-1.31338e-10,1.96852e-14,8201.26,-42.9842], Tmin=(1650.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.5659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=COJ) + radical(C2CsJO)"""),
)

species(
    label = 'C=C([O])CC1[CH]O1(9837)',
    structure = SMILES('C=C([O])CC1[CH]O1'),
    E0 = (83.7599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.310454,0.0705069,-7.04821e-05,3.54012e-08,-6.50764e-12,10250.1,28.958], Tmin=(100,'K'), Tmax=(1633.34,'K')), NASAPolynomial(coeffs=[16.1834,0.011664,5.00841e-07,-4.87085e-10,4.3978e-14,7323.17,-51.1739], Tmin=(1633.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(CCsJO)"""),
)

species(
    label = 'C=C1C[CH]C([O])O1(14394)',
    structure = SMILES('C=C1C[CH]C([O])O1'),
    E0 = (75.6229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14796,0.0179784,8.10856e-05,-1.2221e-07,4.88941e-11,9182.73,23.5872], Tmin=(100,'K'), Tmax=(952.185,'K')), NASAPolynomial(coeffs=[16.7193,0.0143359,-3.86723e-06,7.66216e-10,-6.51395e-14,3798.01,-59.7014], Tmin=(952.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.6229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=C1C[CH][CH]OC1(14261)',
    structure = SMILES('O=C1C[CH][CH]OC1'),
    E0 = (38.8692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50787,0.0388792,1.92844e-05,-5.18267e-08,2.13225e-11,4778.44,20.8051], Tmin=(100,'K'), Tmax=(1035.42,'K')), NASAPolynomial(coeffs=[15.7381,0.0214928,-9.98086e-06,2.07623e-09,-1.57523e-13,-183.266,-58.0724], Tmin=(1035.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclohexanone) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'C=C(O)CC=C=O(14475)',
    structure = SMILES('C=C(O)CC=C=O'),
    E0 = (-212.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563763,0.0691003,-7.25253e-05,3.84619e-08,-7.94608e-12,-25403.7,24.4769], Tmin=(100,'K'), Tmax=(1189.8,'K')), NASAPolynomial(coeffs=[16.3619,0.0159878,-5.56487e-06,9.42366e-10,-6.24177e-14,-29163,-54.4854], Tmin=(1189.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C([O])[CH]C=C[O](14476)',
    structure = SMILES('[CH2]C([O])[CH]C=C[O]'),
    E0 = (278.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,261.077,261.115,261.199,261.366],'cm^-1')),
        HinderedRotor(inertia=(0.670483,'amu*angstrom^2'), symmetry=1, barrier=(32.4642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.670332,'amu*angstrom^2'), symmetry=1, barrier=(32.4648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671162,'amu*angstrom^2'), symmetry=1, barrier=(32.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353507,0.067649,-4.40973e-05,-3.54624e-09,9.50929e-12,33670.3,26.6788], Tmin=(100,'K'), Tmax=(957.005,'K')), NASAPolynomial(coeffs=[20.7377,0.0121256,-3.58495e-06,6.35026e-10,-4.76395e-14,28409.7,-77.8674], Tmin=(957.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])[CH]CC=O(10748)',
    structure = SMILES('[CH2][C]([O])[CH]CC=O'),
    E0 = (395.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,308.165,1383.25,1383.25],'cm^-1')),
        HinderedRotor(inertia=(0.00177525,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116014,'amu*angstrom^2'), symmetry=1, barrier=(7.81791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116012,'amu*angstrom^2'), symmetry=1, barrier=(7.81791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116015,'amu*angstrom^2'), symmetry=1, barrier=(7.81796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80675,0.0762053,-0.000111686,9.37678e-08,-3.1662e-11,47631.9,30.0632], Tmin=(100,'K'), Tmax=(815.37,'K')), NASAPolynomial(coeffs=[8.40302,0.0311962,-1.46391e-05,2.77215e-09,-1.90642e-13,46650.6,-3.45541], Tmin=(815.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CJCO) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C[C]=C[O](14477)',
    structure = SMILES('[CH2]C([O])C[C]=C[O]'),
    E0 = (399.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,430.496,430.584,430.608,430.639],'cm^-1')),
        HinderedRotor(inertia=(0.123564,'amu*angstrom^2'), symmetry=1, barrier=(16.2584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123538,'amu*angstrom^2'), symmetry=1, barrier=(16.2577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123587,'amu*angstrom^2'), symmetry=1, barrier=(16.2589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519967,0.0685223,-6.19679e-05,2.3012e-08,-1.57304e-12,48203.9,29.3555], Tmin=(100,'K'), Tmax=(988.473,'K')), NASAPolynomial(coeffs=[17.519,0.0154745,-5.35534e-06,9.40601e-10,-6.54123e-14,44074.3,-56.3477], Tmin=(988.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C[C]([O])C[C]=C[O](14478)',
    structure = SMILES('C[C]([O])C[C]=C[O]'),
    E0 = (364.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,386.15,386.15,386.151,386.151],'cm^-1')),
        HinderedRotor(inertia=(0.130593,'amu*angstrom^2'), symmetry=1, barrier=(13.8187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130595,'amu*angstrom^2'), symmetry=1, barrier=(13.8187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130595,'amu*angstrom^2'), symmetry=1, barrier=(13.8187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600911,0.0710706,-7.56689e-05,4.1202e-08,-8.83226e-12,43992.2,27.8569], Tmin=(100,'K'), Tmax=(1139.88,'K')), NASAPolynomial(coeffs=[15.3564,0.0192916,-7.53157e-06,1.35151e-09,-9.21978e-14,40628.3,-45.2615], Tmin=(1139.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C([O])CC[CH][O](14479)',
    structure = SMILES('[CH]=C([O])CC[CH][O]'),
    E0 = (376.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,322.014,322.028,322.031,2701.74],'cm^-1')),
        HinderedRotor(inertia=(0.00162558,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130329,'amu*angstrom^2'), symmetry=1, barrier=(9.5906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130316,'amu*angstrom^2'), symmetry=1, barrier=(9.59062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.698414,0.0786659,-0.000117316,9.79708e-08,-3.24629e-11,45361.4,29.2488], Tmin=(100,'K'), Tmax=(842.621,'K')), NASAPolynomial(coeffs=[8.96139,0.0299282,-1.36216e-05,2.53153e-09,-1.71697e-13,44306.6,-7.1965], Tmin=(842.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C[CH][C]=O(10747)',
    structure = SMILES('[CH2]C([O])C[CH][C]=O'),
    E0 = (346.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1855,455,950,3025,407.5,1350,352.5,290.835,290.847,3269.61],'cm^-1')),
        HinderedRotor(inertia=(0.347575,'amu*angstrom^2'), symmetry=1, barrier=(20.8673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00199269,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00199306,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24257,'amu*angstrom^2'), symmetry=1, barrier=(74.5734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927044,0.0683989,-7.7314e-05,4.74186e-08,-1.17062e-11,41734.5,28.2945], Tmin=(100,'K'), Tmax=(984.257,'K')), NASAPolynomial(coeffs=[11.797,0.0242241,-9.99241e-06,1.82005e-09,-1.24287e-13,39594.7,-23.9741], Tmin=(984.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CC(C)OJ) + radical(CCJCHO) + radical(CJCO)"""),
)

species(
    label = '[CH]=C([O])C[CH]C[O](14480)',
    structure = SMILES('[CH]=C([O])C[CH]C[O]'),
    E0 = (395.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,180,180,725.221,2514.24],'cm^-1')),
        HinderedRotor(inertia=(0.000526114,'amu*angstrom^2'), symmetry=1, barrier=(2.36004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614208,'amu*angstrom^2'), symmetry=1, barrier=(14.1218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0378634,'amu*angstrom^2'), symmetry=1, barrier=(14.124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28293,0.0631412,-7.12497e-05,4.69436e-08,-1.29133e-11,47701,30.108], Tmin=(100,'K'), Tmax=(872.886,'K')), NASAPolynomial(coeffs=[8.68141,0.0292383,-1.2991e-05,2.4494e-09,-1.70207e-13,46409.4,-4.5797], Tmin=(872.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = 'C[C]([O])C[CH][C]=O(10752)',
    structure = SMILES('C[C]([O])C[CH][C]=O'),
    E0 = (311.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,360,370,350,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,339.345,2498.06,2498.13],'cm^-1')),
        HinderedRotor(inertia=(0.156733,'amu*angstrom^2'), symmetry=1, barrier=(12.8076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156773,'amu*angstrom^2'), symmetry=1, barrier=(12.8074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156756,'amu*angstrom^2'), symmetry=1, barrier=(12.8074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835902,'amu*angstrom^2'), symmetry=1, barrier=(68.3006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870113,0.0726275,-9.72116e-05,7.42599e-08,-2.30041e-11,37528.8,27.2875], Tmin=(100,'K'), Tmax=(815.035,'K')), NASAPolynomial(coeffs=[9.61238,0.028076,-1.21883e-05,2.23561e-09,-1.51463e-13,36158.4,-12.7654], Tmin=(815.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C1CC=COO1(14351)',
    structure = SMILES('C=C1CC=COO1'),
    E0 = (62.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25607,0.0160686,8.64571e-05,-1.19331e-07,4.45326e-11,7616.99,18.0563], Tmin=(100,'K'), Tmax=(996.899,'K')), NASAPolynomial(coeffs=[14.3379,0.0239958,-1.0341e-05,2.15832e-09,-1.67321e-13,2405.34,-54.2514], Tmin=(996.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane)"""),
)

species(
    label = '[O][C]1C[CH]C([O])C1(14481)',
    structure = SMILES('[O][C]1C[CH]C([O])C1'),
    E0 = (400.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79383,0.0329907,3.16654e-05,-6.73243e-08,2.88597e-11,48279,25.7797], Tmin=(100,'K'), Tmax=(955.784,'K')), NASAPolynomial(coeffs=[14.5393,0.0182526,-5.78665e-06,1.05517e-09,-7.88292e-14,44079.4,-44.3573], Tmin=(955.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJCO)"""),
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
    label = 'C=[C]CC=C[O](13793)',
    structure = SMILES('C=[C]CC=C[O]'),
    E0 = (259.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,193.592,193.612,193.649],'cm^-1')),
        HinderedRotor(inertia=(0.750706,'amu*angstrom^2'), symmetry=1, barrier=(19.9642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750429,'amu*angstrom^2'), symmetry=1, barrier=(19.9636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58152,0.0414199,-1.96358e-06,-3.23802e-08,1.68456e-11,31279.7,23.247], Tmin=(100,'K'), Tmax=(955.649,'K')), NASAPolynomial(coeffs=[14.9255,0.0142834,-4.44361e-06,7.93615e-10,-5.85182e-14,27418,-47.3852], Tmin=(955.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC(=C)[O](9696)',
    structure = SMILES('[CH]=CCC(=C)[O]'),
    E0 = (259.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3010,987.5,1337.5,450,1655,356.919,357.037],'cm^-1')),
        HinderedRotor(inertia=(0.12101,'amu*angstrom^2'), symmetry=1, barrier=(10.9442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121048,'amu*angstrom^2'), symmetry=1, barrier=(10.9435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51949,0.0464528,-2.33308e-05,-6.76514e-09,7.1228e-12,31258.2,23.688], Tmin=(100,'K'), Tmax=(965.886,'K')), NASAPolynomial(coeffs=[13.3363,0.0162286,-5.45376e-06,9.53765e-10,-6.66555e-14,28102.5,-37.4299], Tmin=(965.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O][CH]C1CC(=O)C1(14482)',
    structure = SMILES('[O][CH]C1CC(=O)C1'),
    E0 = (97.4243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03193,0.0438309,-2.28874e-05,5.30365e-09,-4.72265e-13,11787.6,22.9216], Tmin=(100,'K'), Tmax=(2599.32,'K')), NASAPolynomial(coeffs=[20.5984,0.0152609,-6.40108e-06,1.07547e-09,-6.56206e-14,2135.09,-84.3873], Tmin=(2599.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.4243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CC([O])=CC=C[O](14483)',
    structure = SMILES('CC([O])=CC=C[O]'),
    E0 = (-86.1882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756235,0.0562742,-1.51561e-05,-3.6247e-08,2.26938e-11,-10235.2,24.2724], Tmin=(100,'K'), Tmax=(917.521,'K')), NASAPolynomial(coeffs=[21.0957,0.00729401,3.01004e-08,-1.33448e-10,7.27657e-15,-15638.3,-81.2071], Tmin=(917.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.1882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)C[C]=C[O](14484)',
    structure = SMILES('CC(=O)C[C]=C[O]'),
    E0 = (22.5022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,480.367,482.586,483.3],'cm^-1')),
        HinderedRotor(inertia=(0.100059,'amu*angstrom^2'), symmetry=1, barrier=(16.3483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100012,'amu*angstrom^2'), symmetry=1, barrier=(16.3359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996738,'amu*angstrom^2'), symmetry=1, barrier=(16.3571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16409,0.0509271,-2.00839e-05,-1.09955e-08,7.32876e-12,2818.23,25.1855], Tmin=(100,'K'), Tmax=(1083.49,'K')), NASAPolynomial(coeffs=[15.6593,0.0202904,-9.34022e-06,1.88049e-09,-1.38423e-13,-1665.63,-52.1043], Tmin=(1083.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.5022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C[CH][C]=O(14485)',
    structure = SMILES('CC(=O)C[CH][C]=O'),
    E0 = (-4.91091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,236.127,3364.26],'cm^-1')),
        HinderedRotor(inertia=(0.131212,'amu*angstrom^2'), symmetry=1, barrier=(5.18918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131221,'amu*angstrom^2'), symmetry=1, barrier=(5.18968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00302493,'amu*angstrom^2'), symmetry=1, barrier=(0.11963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131179,'amu*angstrom^2'), symmetry=1, barrier=(5.18849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31242,0.0472865,6.26782e-06,-1.42716e-07,1.51291e-10,-540.757,24.8893], Tmin=(100,'K'), Tmax=(424.196,'K')), NASAPolynomial(coeffs=[5.91088,0.033401,-1.55181e-05,2.92749e-09,-2.00695e-13,-1026.41,8.48874], Tmin=(424.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.91091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCC=O)"""),
)

species(
    label = '[O]C1[CH]CC(=O)C1(14305)',
    structure = SMILES('[O]C1[CH]CC(=O)C1'),
    E0 = (43.0727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0513,0.0278405,4.08901e-05,-7.07999e-08,2.83012e-11,5264.05,22.6876], Tmin=(100,'K'), Tmax=(977.732,'K')), NASAPolynomial(coeffs=[12.2898,0.0229881,-8.48256e-06,1.6055e-09,-1.18159e-13,1491.76,-35.5293], Tmin=(977.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C1C[CH][CH]OO1(14342)',
    structure = SMILES('C=C1C[CH][CH]OO1'),
    E0 = (321.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77053,0.0338397,2.89147e-05,-5.56262e-08,2.08717e-11,38749.5,21.9801], Tmin=(100,'K'), Tmax=(1074.43,'K')), NASAPolynomial(coeffs=[13.6238,0.0268861,-1.32774e-05,2.75659e-09,-2.06017e-13,34056.7,-46.0415], Tmin=(1074.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(CCJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = 'CC(=O)CC=C=O(14486)',
    structure = SMILES('CC(=O)CC=C=O'),
    E0 = (-237.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82806,0.048822,-3.3564e-05,1.1882e-08,-1.77819e-12,-28495.5,23.1054], Tmin=(100,'K'), Tmax=(1476.21,'K')), NASAPolynomial(coeffs=[9.50833,0.0280114,-1.24183e-05,2.33256e-09,-1.60977e-13,-30763,-16.9388], Tmin=(1476.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = '[CH2][C](O)C[C]=C[O](14487)',
    structure = SMILES('[CH2][C](O)C[C]=C[O]'),
    E0 = (345.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,297.394,297.407,297.41],'cm^-1')),
        HinderedRotor(inertia=(0.00190587,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226323,'amu*angstrom^2'), symmetry=1, barrier=(14.2059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226319,'amu*angstrom^2'), symmetry=1, barrier=(14.2058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22633,'amu*angstrom^2'), symmetry=1, barrier=(14.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341857,0.0795728,-9.88872e-05,6.22261e-08,-1.5265e-11,41741.5,29.6142], Tmin=(100,'K'), Tmax=(1004.23,'K')), NASAPolynomial(coeffs=[15.7824,0.0180715,-7.02493e-06,1.2433e-09,-8.36837e-14,38640.3,-44.943], Tmin=(1004.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(Cds_S) + radical(C=COJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C](O)C[CH][C]=O(10751)',
    structure = SMILES('[CH2][C](O)C[CH][C]=O'),
    E0 = (292.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1277.5,1000,3025,407.5,1350,352.5,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,342.756,3087.5],'cm^-1')),
        HinderedRotor(inertia=(0.105524,'amu*angstrom^2'), symmetry=1, barrier=(8.79732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105524,'amu*angstrom^2'), symmetry=1, barrier=(8.79732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143492,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105524,'amu*angstrom^2'), symmetry=1, barrier=(8.79732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812491,'amu*angstrom^2'), symmetry=1, barrier=(67.7362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497798,0.0825178,-0.000125476,1.01984e-07,-3.23111e-11,35282.8,29.447], Tmin=(100,'K'), Tmax=(881,'K')), NASAPolynomial(coeffs=[10.7297,0.0256701,-1.09931e-05,1.96389e-09,-1.2933e-13,33883.2,-16.3308], Tmin=(881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CCJCHO) + radical(C2CsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CCC[C]=O(12768)',
    structure = SMILES('[O]C=CCC[C]=O'),
    E0 = (-30.4708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,323.809,324.268,324.608],'cm^-1')),
        HinderedRotor(inertia=(0.180073,'amu*angstrom^2'), symmetry=1, barrier=(13.5281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179602,'amu*angstrom^2'), symmetry=1, barrier=(13.5261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181535,'amu*angstrom^2'), symmetry=1, barrier=(13.528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4066.4,'J/mol'), sigma=(6.53709,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.16 K, Pc=33.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.946714,0.0595963,-4.7233e-05,1.46353e-08,-1.92794e-13,-3548.27,27.4017], Tmin=(100,'K'), Tmax=(1020.5,'K')), NASAPolynomial(coeffs=[14.848,0.0186152,-6.85012e-06,1.22432e-09,-8.47827e-14,-7088.86,-43.3922], Tmin=(1020.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.4708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC=COC1(14255)',
    structure = SMILES('O=C1CC=COC1'),
    E0 = (-263.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01821,0.0107983,0.000126773,-1.75853e-07,6.73814e-11,-31595.1,13.9998], Tmin=(100,'K'), Tmax=(978.226,'K')), NASAPolynomial(coeffs=[21.5044,0.0142334,-5.94096e-06,1.44884e-09,-1.27841e-13,-39384.2,-99.9071], Tmin=(978.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclohexane)"""),
)

species(
    label = '[CH2]C(=C)OC=C[O](14488)',
    structure = SMILES('[CH2]C(=C)OC=C[O]'),
    E0 = (19.6089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,310.07,310.132,310.338,310.339],'cm^-1')),
        HinderedRotor(inertia=(0.322234,'amu*angstrom^2'), symmetry=1, barrier=(22.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322962,'amu*angstrom^2'), symmetry=1, barrier=(22.0317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322751,'amu*angstrom^2'), symmetry=1, barrier=(22.0308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722334,0.0577635,-2.01052e-05,-2.95389e-08,1.97202e-11,2489.68,26.0745], Tmin=(100,'K'), Tmax=(924.909,'K')), NASAPolynomial(coeffs=[20.5909,0.00892866,-1.06013e-06,9.23534e-11,-8.71335e-15,-2772.14,-76.8053], Tmin=(924.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.6089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][C]1C[CH]C([O])O1(14489)',
    structure = SMILES('[CH2][C]1C[CH]C([O])O1'),
    E0 = (381.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87442,0.0478901,-2.94147e-05,5.05538e-10,7.15437e-12,45972.2,23.285], Tmin=(100,'K'), Tmax=(707.263,'K')), NASAPolynomial(coeffs=[6.92334,0.0290421,-1.00272e-05,1.63558e-09,-1.04142e-13,45015.3,-1.04062], Tmin=(707.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(C2CsJOCs) + radical(CCOJ) + radical(CJC(C)OC)"""),
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
    label = '[O]C=CC[C]=O(4679)',
    structure = SMILES('[O]C=CC[C]=O'),
    E0 = (-0.508181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,189.447,189.447],'cm^-1')),
        HinderedRotor(inertia=(0.845448,'amu*angstrom^2'), symmetry=1, barrier=(21.5313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845432,'amu*angstrom^2'), symmetry=1, barrier=(21.5313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88818,0.0328175,9.85639e-06,-4.13901e-08,1.88828e-11,27.0453,21.3343], Tmin=(100,'K'), Tmax=(992.314,'K')), NASAPolynomial(coeffs=[15.7297,0.00949103,-3.96267e-06,8.67221e-10,-7.04814e-14,-4318.53,-53.3907], Tmin=(992.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.508181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1([O])CC1C=O(14397)',
    structure = SMILES('[CH2]C1([O])CC1C=O'),
    E0 = (156.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796336,0.0624066,-5.54408e-05,2.50035e-08,-4.46098e-12,18936.3,24.1339], Tmin=(100,'K'), Tmax=(1356.2,'K')), NASAPolynomial(coeffs=[15.5991,0.0187475,-7.15305e-06,1.26694e-09,-8.54646e-14,14921.1,-51.7913], Tmin=(1356.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C=CC=O(14490)',
    structure = SMILES('C=C([O])C=CC=O'),
    E0 = (-105.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.804431,'amu*angstrom^2'), symmetry=1, barrier=(18.4954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804525,'amu*angstrom^2'), symmetry=1, barrier=(18.4976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16866,0.0608176,-5.94524e-05,2.9485e-08,-5.85639e-12,-12604.2,21.9805], Tmin=(100,'K'), Tmax=(1208.63,'K')), NASAPolynomial(coeffs=[13.3958,0.0203512,-9.23039e-06,1.78305e-09,-1.26336e-13,-15559.8,-39.3255], Tmin=(1208.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([O])=CCC=O(10755)',
    structure = SMILES('[CH2]C([O])=CCC=O'),
    E0 = (-47.0015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,464.496,465.158],'cm^-1')),
        HinderedRotor(inertia=(0.000782822,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947177,'amu*angstrom^2'), symmetry=1, barrier=(14.4834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0949332,'amu*angstrom^2'), symmetry=1, barrier=(14.486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16183,0.0531123,-3.14972e-05,2.75316e-09,2.39166e-12,-5542.96,25.9129], Tmin=(100,'K'), Tmax=(1122.95,'K')), NASAPolynomial(coeffs=[14.5717,0.0208532,-9.121e-06,1.76656e-09,-1.26489e-13,-9532.44,-44.6904], Tmin=(1122.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.0015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])CC[C]=O(9776)',
    structure = SMILES('C=C([O])CC[C]=O'),
    E0 = (-39.8946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,313.595,313.599,4000],'cm^-1')),
        HinderedRotor(inertia=(0.176709,'amu*angstrom^2'), symmetry=1, barrier=(12.3299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176675,'amu*angstrom^2'), symmetry=1, barrier=(12.3299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176672,'amu*angstrom^2'), symmetry=1, barrier=(12.3299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15879,0.0628988,-6.6679e-05,3.88165e-08,-9.17427e-12,-4696.2,26.9339], Tmin=(100,'K'), Tmax=(1022.11,'K')), NASAPolynomial(coeffs=[11.036,0.0242448,-9.95222e-06,1.8167e-09,-1.24421e-13,-6715.32,-20.9338], Tmin=(1022.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.8946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])CCC=O(14491)',
    structure = SMILES('[CH]=C([O])CCC=O'),
    E0 = (47.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2782.5,750,1395,475,1775,1000,199.962,200.141],'cm^-1')),
        HinderedRotor(inertia=(0.00419835,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222636,'amu*angstrom^2'), symmetry=1, barrier=(6.35104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223353,'amu*angstrom^2'), symmetry=1, barrier=(6.35336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22681,0.0626074,-6.5164e-05,3.76681e-08,-8.96047e-12,5780.32,26.2769], Tmin=(100,'K'), Tmax=(1007.86,'K')), NASAPolynomial(coeffs=[10.3177,0.0265267,-1.14638e-05,2.1463e-09,-1.49082e-13,3947.89,-17.6523], Tmin=(1007.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C=CC=O(14492)',
    structure = SMILES('C=C(O)C=CC=O'),
    E0 = (-243.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.402004,0.0694042,-6.75524e-05,3.19153e-08,-5.84626e-12,-29143.7,23.2376], Tmin=(100,'K'), Tmax=(1337.05,'K')), NASAPolynomial(coeffs=[18.8569,0.0141932,-5.61246e-06,1.03134e-09,-7.16018e-14,-34078.7,-71.1572], Tmin=(1337.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[CH2]C(C=O)C(=C)[O](12762)',
    structure = SMILES('[CH2]C(C=O)C(=C)[O]'),
    E0 = (3.76026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,290.147,290.154],'cm^-1')),
        HinderedRotor(inertia=(0.110841,'amu*angstrom^2'), symmetry=1, barrier=(6.62188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110844,'amu*angstrom^2'), symmetry=1, barrier=(6.62184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110844,'amu*angstrom^2'), symmetry=1, barrier=(6.62203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4030.31,'J/mol'), sigma=(6.52518,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.53 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.912788,0.071282,-9.01776e-05,6.36277e-08,-1.81958e-11,560.377,26.5475], Tmin=(100,'K'), Tmax=(851.271,'K')), NASAPolynomial(coeffs=[10.3893,0.0267531,-1.17142e-05,2.17942e-09,-1.4976e-13,-1053.04,-17.6453], Tmin=(851.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.76026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1CC(C=O)O1(14327)',
    structure = SMILES('C=C1CC(C=O)O1'),
    E0 = (-171.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7373,0.0243213,7.89174e-05,-1.35208e-07,5.84286e-11,-20566.1,19.425], Tmin=(100,'K'), Tmax=(914.759,'K')), NASAPolynomial(coeffs=[21.4272,0.00443927,2.93902e-06,-7.02941e-10,4.24279e-14,-26938.9,-88.956], Tmin=(914.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
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
    label = '[CH]CC(=C)[O](4317)',
    structure = SMILES('[CH]CC(=C)[O]'),
    E0 = (354.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180,815.109,819.63,819.773,820.774],'cm^-1')),
        HinderedRotor(inertia=(0.00156324,'amu*angstrom^2'), symmetry=1, barrier=(14.8523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647586,'amu*angstrom^2'), symmetry=1, barrier=(14.8893,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96185,0.0381142,-1.88956e-05,-6.94752e-09,6.76872e-12,42768.9,19.9557], Tmin=(100,'K'), Tmax=(946.625,'K')), NASAPolynomial(coeffs=[11.9364,0.0122826,-3.81816e-06,6.42341e-10,-4.4461e-14,40149.4,-31.4801], Tmin=(946.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC(=O)C=CC=O(14493)',
    structure = SMILES('CC(=O)C=CC=O'),
    E0 = (-249.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72359,0.0430426,-2.26393e-05,4.45245e-09,-2.82809e-13,-30034.5,19.0288], Tmin=(100,'K'), Tmax=(2443.83,'K')), NASAPolynomial(coeffs=[29.6812,0.00674437,-5.16283e-06,9.95166e-10,-6.31701e-14,-45547.2,-139.895], Tmin=(2443.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=CC1CC(=O)C1(12774)',
    structure = SMILES('O=CC1CC(=O)C1'),
    E0 = (-228.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04545,0.034498,6.21942e-06,-2.48178e-08,9.6488e-12,-27377,21.925], Tmin=(100,'K'), Tmax=(1096.75,'K')), NASAPolynomial(coeffs=[8.83253,0.0288784,-1.22637e-05,2.32414e-09,-1.6399e-13,-30016.5,-16.6918], Tmin=(1096.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Cyclobutanone)"""),
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
    E0 = (-55.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (154.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (65.3864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (140.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (76.7207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (205.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (114.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (246.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (191.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (89.0115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (59.0984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (34.1279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (168.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (251.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (310.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (273.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (394.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (404.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (532.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (175.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (158.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (75.6229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (57.6842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-30.3608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (301.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (417.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (488.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (403.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (401.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (371.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (420.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (336.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (62.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (550.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (666.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (665.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (97.4243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (64.8625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (66.8108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (34.9884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (43.0727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (321.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-30.3608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (385.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (317.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (215.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-47.1334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (333.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (531.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (415.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (157.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (115.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (84.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (69.0669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (116.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (91.5495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (22.9131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (163.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (-47.4263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (422.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (22.9131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (-47.4263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C=O(598)', 'C=CC=O(5269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C1CC([CH][O])O1(14462)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(209.949,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['[CH2]C1([O])CC=CO1(14440)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(120.72,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 119.5 to 120.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C=C([O])CC=C=O(14463)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=O(598)', '[CH2]C=C[O](5266)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000178913,'m^3/(mol*s)'), n=2.74787, Ea=(45.8259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CsJ-CdHH] for rate rule [Ck_O;CsJ-CdHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])CC=[C]O(14464)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C(O)[CH]C=C[O](14465)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(O)CC=C[O](14466)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C[C]=CO(14467)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C[C]=C[O](14468)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C([O])[CH]C=CO(14469)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C[CH][C]=O(14470)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])CC=CO(14471)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(173703,'s^-1'), n=1.89007, Ea=(118.15,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(601)', '[CH2]C=C[O](5266)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C[O](602)', '[CH2]C(=C)[O](4273)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.33799e+08,'m^3/(mol*s)'), n=-0.455312, Ea=(0.377199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cd] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', 'C=C([O])[CH]C=C[O](14472)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.14e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', 'C=C([O])C[C]=C[O](14473)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]=C([O])CC=C[O](14474)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2][C]([O])CC=C=O(10744)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['[O]C=CC[C]1CO1(9768)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C([O])CC1[CH]O1(9837)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C1C[CH]C([O])O1(14394)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(130.957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 126.4 to 131.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['O=C1C[CH][CH]OC1(14261)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(113.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C(O)CC=C=O(14475)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([O])[CH]C=C[O](14476)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]([O])[CH]CC=O(10748)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([O])C[C]=C[O](14477)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[C]([O])C[C]=C[O](14478)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C([O])CC[CH][O](14479)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([O])C[CH][C]=O(10747)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([O])C[CH]C[O](14480)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[C]([O])C[CH][C]=O(10752)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C1CC=COO1(14351)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(117.979,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination
Ea raised from 113.4 to 118.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O][C]1C[CH]C([O])C1(14481)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', 'C=[C]CC=C[O](13793)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH]=CCC(=C)[O](9696)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['[O][CH]C1CC(=O)C1(14482)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(177207,'s^-1'), n=1.88643, Ea=(152.758,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 151.6 to 152.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['CC([O])=CC=C[O](14483)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.00568695,'s^-1'), n=4.30267, Ea=(120.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC(=O)C[C]=C[O](14484)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC(=O)C[CH][C]=O(14485)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['[O]C1[CH]CC(=O)C1(14305)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.47079e+07,'s^-1'), n=0.909323, Ea=(98.4068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 95.0 to 98.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C1C[CH][CH]OO1(14342)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.39072e+10,'s^-1'), n=0.346137, Ea=(376.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Endocyclic
Ea raised from 374.3 to 376.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['CC(=O)CC=C=O(14486)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][C](O)C[C]=C[O](14487)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C](O)C[CH][C]=O(10751)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=CCC[C]=O(12768)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['O=C1CC=COC1(14255)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Ypri_rad_out] + [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(=C)OC=C[O](14488)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][C]1C[CH]C([O])O1(14489)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction50',
    reactants = ['CH2(T)(28)', '[O]C=CC[C]=O(4679)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['[CH2]C1([O])CC1C=O(14397)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(212.728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['H(8)', 'C=C([O])C=CC=O(14490)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(30.9498,'m^3/(mol*s)'), n=1.77449, Ea=(9.26302,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-OneDeH;HJ] for rate rule [Cds-CdH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][C]=O(601)', 'C=CC=O(5269)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(0.0102751,'m^3/(mol*s)'), n=2.40501, Ea=(4.48561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CJ] for rate rule [Cds-HH_Cds-COH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C([O])=CCC=O(10755)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C([O])CC[C]=O(9776)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=C([O])CCC=O(14491)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C(O)C=CC=O(14492)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=O)C(=C)[O](12762)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['C=C1CC(C=O)O1(14327)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]=O(373)', '[CH]CC(=C)[O](4317)'],
    products = ['C=C([O])CC=C[O](12760)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['CC(=O)C=CC=O(14493)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C([O])CC=C[O](12760)'],
    products = ['O=CC1CC(=O)C1(12774)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '3519',
    isomers = [
        'C=C([O])CC=C[O](12760)',
    ],
    reactants = [
        ('C=C=O(598)', 'C=CC=O(5269)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3519',
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


species(
    label = '[CH2]C(=CO)C[C]=O(14236)',
    structure = SMILES('[CH2]C(=CO)C[C]=O'),
    E0 = (-29.5266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.954771,'amu*angstrom^2'), symmetry=1, barrier=(21.9521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952009,'amu*angstrom^2'), symmetry=1, barrier=(21.8886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951705,'amu*angstrom^2'), symmetry=1, barrier=(21.8816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.953727,'amu*angstrom^2'), symmetry=1, barrier=(21.9281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636097,0.057278,-1.55871e-05,-3.36923e-08,2.00748e-11,-3415,25.454], Tmin=(100,'K'), Tmax=(963.153,'K')), NASAPolynomial(coeffs=[22.0242,0.00863463,-2.40919e-06,5.0148e-10,-4.37789e-14,-9398.75,-86.6034], Tmin=(963.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Allyl_P)"""),
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
    label = '[CH2]C1([CH]O)CC1=O(28256)',
    structure = SMILES('[CH2]C1([CH]O)CC1=O'),
    E0 = (166.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303633,0.0816349,-0.000106952,7.22756e-08,-1.90679e-11,20202.3,24.5054], Tmin=(100,'K'), Tmax=(935.309,'K')), NASAPolynomial(coeffs=[14.8338,0.019494,-7.29247e-06,1.24027e-09,-8.06009e-14,17484.3,-44.6223], Tmin=(935.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CJC(C)2C=O) + radical(CCsJOH)"""),
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
    label = '[CH2]C(C=C=O)=CO(28257)',
    structure = SMILES('[CH2]C(C=C=O)=CO'),
    E0 = (-64.3733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.13456,'amu*angstrom^2'), symmetry=1, barrier=(26.0857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13549,'amu*angstrom^2'), symmetry=1, barrier=(26.1071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13377,'amu*angstrom^2'), symmetry=1, barrier=(26.0675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0297091,0.0727639,-7.83263e-05,4.0091e-08,-7.65739e-12,-7583.65,26.4283], Tmin=(100,'K'), Tmax=(1459.78,'K')), NASAPolynomial(coeffs=[20.5436,0.00662746,-3.3602e-07,-1.07812e-10,1.16113e-14,-12550,-77.0458], Tmin=(1459.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.3733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH) + radical(C=C(CJ)C=C=O)"""),
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
    label = '[CH2]C(C=C[O])=CO(28258)',
    structure = SMILES('[CH2]C(C=C[O])=CO'),
    E0 = (-68.2929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.181141,0.0588739,1.04657e-05,-8.65102e-08,4.66728e-11,-8052.6,23.7149], Tmin=(100,'K'), Tmax=(902.471,'K')), NASAPolynomial(coeffs=[30.5636,-0.00790706,8.63497e-06,-1.81064e-09,1.2101e-13,-16300.8,-135.061], Tmin=(902.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.2929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC([CH]O)=C[C]=O(28259)',
    structure = SMILES('CC([CH]O)=C[C]=O'),
    E0 = (-31.6386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15985,0.0657526,-6.63459e-05,3.52034e-08,-7.69944e-12,-3705.68,25.0797], Tmin=(100,'K'), Tmax=(1082.86,'K')), NASAPolynomial(coeffs=[11.4817,0.0276255,-1.35327e-05,2.68959e-09,-1.93151e-13,-5941.15,-25.539], Tmin=(1082.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.6386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC(=[C]O)C[C]=O(28260)',
    structure = SMILES('CC(=[C]O)C[C]=O'),
    E0 = (58.7184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649068,0.0636369,-5.73331e-05,2.50576e-08,-4.26932e-12,7191.25,28.2694], Tmin=(100,'K'), Tmax=(1427.09,'K')), NASAPolynomial(coeffs=[18.0306,0.0149181,-6.12521e-06,1.1358e-09,-7.86665e-14,2230.24,-61.7682], Tmin=(1427.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.7184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O)CC=O(28261)',
    structure = SMILES('[CH2]C(=[C]O)CC=O'),
    E0 = (50.2569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2782.5,750,1395,475,1775,1000,446.416],'cm^-1')),
        HinderedRotor(inertia=(0.143133,'amu*angstrom^2'), symmetry=1, barrier=(19.7799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141102,'amu*angstrom^2'), symmetry=1, barrier=(19.8036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142112,'amu*angstrom^2'), symmetry=1, barrier=(19.809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13987,'amu*angstrom^2'), symmetry=1, barrier=(19.7781,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882609,0.0564477,-2.84763e-05,-8.26108e-09,7.91633e-12,6167.3,26.882], Tmin=(100,'K'), Tmax=(1032.73,'K')), NASAPolynomial(coeffs=[17.5834,0.01677,-7.16993e-06,1.4333e-09,-1.06794e-13,1384.2,-60.6844], Tmin=(1032.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.2569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = 'CC(=C[O])C[C]=O(14220)',
    structure = SMILES('CC(=C[O])C[C]=O'),
    E0 = (-39.5632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,302.456,302.549],'cm^-1')),
        HinderedRotor(inertia=(0.273748,'amu*angstrom^2'), symmetry=1, barrier=(17.7659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27441,'amu*angstrom^2'), symmetry=1, barrier=(17.7669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274448,'amu*angstrom^2'), symmetry=1, barrier=(17.7672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03017,0.0529583,-2.10177e-05,-1.41735e-08,9.5504e-12,-4640.69,25.1724], Tmin=(100,'K'), Tmax=(1035.45,'K')), NASAPolynomial(coeffs=[16.9228,0.0176425,-7.63549e-06,1.53324e-09,-1.14329e-13,-9329.9,-58.8047], Tmin=(1035.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.5632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C[O])CC=O(14224)',
    structure = SMILES('[CH2]C(=C[O])CC=O'),
    E0 = (-48.0246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2782.5,750,1395,475,1775,1000,401.037,401.097],'cm^-1')),
        HinderedRotor(inertia=(0.223666,'amu*angstrom^2'), symmetry=1, barrier=(25.5261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223637,'amu*angstrom^2'), symmetry=1, barrier=(25.5269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223576,'amu*angstrom^2'), symmetry=1, barrier=(25.5265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10933,0.0474297,2.78597e-06,-4.20203e-08,1.98606e-11,-5657.7,24.3495], Tmin=(100,'K'), Tmax=(1006.64,'K')), NASAPolynomial(coeffs=[18.2673,0.0166705,-7.13968e-06,1.48157e-09,-1.14392e-13,-11008,-67.9586], Tmin=(1006.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = '[CH]C(=C)C[C]=O(17871)',
    structure = SMILES('[CH]C(=C)C[C]=O'),
    E0 = (398.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,409.641,409.666,409.679,409.693],'cm^-1')),
        HinderedRotor(inertia=(0.412124,'amu*angstrom^2'), symmetry=1, barrier=(49.0834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412124,'amu*angstrom^2'), symmetry=1, barrier=(49.0829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412112,'amu*angstrom^2'), symmetry=1, barrier=(49.0831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3453.98,'J/mol'), sigma=(5.80287,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.50 K, Pc=40.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6449,0.0440109,-1.8193e-05,-3.18709e-09,2.77889e-12,48013.7,22.8183], Tmin=(100,'K'), Tmax=(1230.58,'K')), NASAPolynomial(coeffs=[11.7847,0.0247584,-1.14335e-05,2.20259e-09,-1.54966e-13,44480.3,-32.4214], Tmin=(1230.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C[O])C[C]=O(14226)',
    structure = SMILES('[CH2]C(=C[O])C[C]=O'),
    E0 = (111.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,370.932,372.05],'cm^-1')),
        HinderedRotor(inertia=(0.236743,'amu*angstrom^2'), symmetry=1, barrier=(23.3386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239156,'amu*angstrom^2'), symmetry=1, barrier=(23.3108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236257,'amu*angstrom^2'), symmetry=1, barrier=(23.3336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05074,0.0502439,-1.18185e-05,-2.81841e-08,1.57501e-11,13582,25.3774], Tmin=(100,'K'), Tmax=(995.541,'K')), NASAPolynomial(coeffs=[18.8748,0.0121531,-4.93878e-06,1.03464e-09,-8.16206e-14,8371.72,-68.8775], Tmin=(995.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]O)=C[C]=O(28262)',
    structure = SMILES('[CH2]C([CH]O)=C[C]=O'),
    E0 = (119.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1855,455,950,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,485.468],'cm^-1')),
        HinderedRotor(inertia=(0.129108,'amu*angstrom^2'), symmetry=1, barrier=(21.5855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93892,'amu*angstrom^2'), symmetry=1, barrier=(21.5876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93891,'amu*angstrom^2'), symmetry=1, barrier=(21.5874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149648,'amu*angstrom^2'), symmetry=1, barrier=(21.5876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959078,0.0656959,-6.65776e-05,3.346e-08,-6.69151e-12,14526.4,26.0745], Tmin=(100,'K'), Tmax=(1203.66,'K')), NASAPolynomial(coeffs=[14.7184,0.0199707,-9.59489e-06,1.89908e-09,-1.36309e-13,11214.1,-42.857], Tmin=(1203.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]O)C[C]=O(28263)',
    structure = SMILES('[CH2]C(=[C]O)C[C]=O'),
    E0 = (210.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,406.632],'cm^-1')),
        HinderedRotor(inertia=(0.155676,'amu*angstrom^2'), symmetry=1, barrier=(18.259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155686,'amu*angstrom^2'), symmetry=1, barrier=(18.2594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(18.2592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155556,'amu*angstrom^2'), symmetry=1, barrier=(18.2608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825085,0.0592546,-4.30832e-05,5.62279e-09,3.76413e-12,25406.9,27.9058], Tmin=(100,'K'), Tmax=(1019.25,'K')), NASAPolynomial(coeffs=[18.1557,0.0123095,-5.00061e-06,9.93631e-10,-7.46112e-14,20779.7,-61.4032], Tmin=(1019.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = 'O=[C]C[C]1CC1O(28264)',
    structure = SMILES('O=[C]C[C]1CC1O'),
    E0 = (42.7353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50993,0.0450551,-1.18e-05,-1.76944e-08,1.01757e-11,5238.1,25.4638], Tmin=(100,'K'), Tmax=(996.434,'K')), NASAPolynomial(coeffs=[12.9948,0.020644,-7.70787e-06,1.41591e-09,-1.00572e-13,1872.42,-35.3067], Tmin=(996.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.7353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCCJ=O) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1CC(=O)C1O(28265)',
    structure = SMILES('[CH2][C]1CC(=O)C1O'),
    E0 = (39.0234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19656,0.0199505,6.9685e-05,-1.07467e-07,4.33726e-11,4776.38,25.9926], Tmin=(100,'K'), Tmax=(945.034,'K')), NASAPolynomial(coeffs=[14.7372,0.0168014,-4.56971e-06,8.24297e-10,-6.48928e-14,176.449,-45.5967], Tmin=(945.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.0234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC(C=C=O)=CO(28266)',
    structure = SMILES('CC(C=C=O)=CO'),
    E0 = (-222.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.382857,0.0736219,-7.76845e-05,3.75369e-08,-5.89298e-12,-26617,23.6883], Tmin=(100,'K'), Tmax=(950.723,'K')), NASAPolynomial(coeffs=[17.7076,0.0139298,-4.3292e-06,7.0079e-10,-4.63632e-14,-30507.7,-62.1549], Tmin=(950.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH)"""),
)

species(
    label = '[CH2][C](C[O])C[C]=O(9762)',
    structure = SMILES('[CH2][C](C[O])C[C]=O'),
    E0 = (323.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,180,1508.52,1508.57],'cm^-1')),
        HinderedRotor(inertia=(0.0852477,'amu*angstrom^2'), symmetry=1, barrier=(1.96001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865277,'amu*angstrom^2'), symmetry=1, barrier=(1.98944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0862898,'amu*angstrom^2'), symmetry=1, barrier=(1.98397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868357,'amu*angstrom^2'), symmetry=1, barrier=(1.99652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34552,0.0664319,-9.86809e-05,9.07187e-08,-3.27591e-11,39004.2,29.488], Tmin=(100,'K'), Tmax=(858.231,'K')), NASAPolynomial(coeffs=[3.72159,0.037553,-1.70887e-05,3.16605e-09,-2.13965e-13,39252.1,22.2082], Tmin=(858.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C([CH]O)[CH][C]=O(14234)',
    structure = SMILES('[CH2]C([CH]O)[CH][C]=O'),
    E0 = (293.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,1855,455,950,1380,1390,370,380,2900,435,243.369,243.369],'cm^-1')),
        HinderedRotor(inertia=(0.00883654,'amu*angstrom^2'), symmetry=1, barrier=(79.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263186,'amu*angstrom^2'), symmetry=1, barrier=(11.0617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00883657,'amu*angstrom^2'), symmetry=1, barrier=(79.5431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284618,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89254,'amu*angstrom^2'), symmetry=1, barrier=(79.5427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.58871,0.0784191,-0.000112216,8.6033e-08,-2.5602e-11,35373.5,30.4583], Tmin=(100,'K'), Tmax=(933.511,'K')), NASAPolynomial(coeffs=[11.2905,0.0232721,-8.67436e-06,1.42744e-09,-8.88874e-14,33780.3,-18.267], Tmin=(933.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCsJOH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=CO)C(=C)[O](14417)',
    structure = SMILES('[CH2]C(=CO)C(=C)[O]'),
    E0 = (-77.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14145,'amu*angstrom^2'), symmetry=1, barrier=(26.2443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14603,'amu*angstrom^2'), symmetry=1, barrier=(26.3495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13543,'amu*angstrom^2'), symmetry=1, barrier=(26.1057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36898,0.0862965,-9.44001e-05,4.72742e-08,-8.53982e-12,-9046.33,28.7692], Tmin=(100,'K'), Tmax=(1645.82,'K')), NASAPolynomial(coeffs=[24.2343,0.00123439,3.93822e-06,-9.89888e-10,7.20786e-14,-14381.2,-98.1124], Tmin=(1645.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=[C]CC[C]=CO(14269)',
    structure = SMILES('O=[C]CC[C]=CO'),
    E0 = (65.9083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,1855,455,950,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.6223,'amu*angstrom^2'), symmetry=1, barrier=(14.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622525,'amu*angstrom^2'), symmetry=1, barrier=(14.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.621562,'amu*angstrom^2'), symmetry=1, barrier=(14.2909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622474,'amu*angstrom^2'), symmetry=1, barrier=(14.3119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4066.4,'J/mol'), sigma=(6.53709,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.16 K, Pc=33.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533613,0.0708145,-7.43059e-05,3.64459e-08,-5.9964e-12,8056.65,27.9472], Tmin=(100,'K'), Tmax=(952.679,'K')), NASAPolynomial(coeffs=[16.6894,0.0147987,-4.71478e-06,7.67731e-10,-5.05565e-14,4442.13,-52.0264], Tmin=(952.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.9083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=C1CC(=CO)C1(28267)',
    structure = SMILES('O=C1CC(=CO)C1'),
    E0 = (-208.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69938,0.0204499,9.53217e-05,-1.46044e-07,5.78953e-11,-24933.1,18.7639], Tmin=(100,'K'), Tmax=(976.416,'K')), NASAPolynomial(coeffs=[22.9466,0.00915151,-3.68072e-06,9.98838e-10,-9.45577e-14,-32692.9,-101.724], Tmin=(976.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + ring(Cyclobutane)"""),
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
    label = '[CH2]C([CH2])=CO(17684)',
    structure = SMILES('[CH2]C([CH2])=CO'),
    E0 = (61.2744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.31967,'amu*angstrom^2'), symmetry=1, barrier=(30.3418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32161,'amu*angstrom^2'), symmetry=1, barrier=(30.3865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32215,'amu*angstrom^2'), symmetry=1, barrier=(30.3989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44673,0.0388097,1.68219e-05,-6.63718e-08,3.33598e-11,7477.77,17.5086], Tmin=(100,'K'), Tmax=(909.843,'K')), NASAPolynomial(coeffs=[20.4722,0.00116961,3.03528e-06,-6.98811e-10,4.52843e-14,2111.65,-82.9443], Tmin=(909.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.2744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=O)C[C]=O(12770)',
    structure = SMILES('[CH2]C(C=O)C[C]=O'),
    E0 = (31.0278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.24716,'amu*angstrom^2'), symmetry=1, barrier=(5.68269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247079,'amu*angstrom^2'), symmetry=1, barrier=(5.68084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247163,'amu*angstrom^2'), symmetry=1, barrier=(5.68276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246997,'amu*angstrom^2'), symmetry=1, barrier=(5.67894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.89,'J/mol'), sigma=(6.24759,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.19 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.886528,0.0762956,-0.000118929,1.05619e-07,-3.68565e-11,3836.37,28.1193], Tmin=(100,'K'), Tmax=(844.531,'K')), NASAPolynomial(coeffs=[6.78841,0.0332514,-1.56737e-05,2.95128e-09,-2.01291e-13,3377.66,3.82956], Tmin=(844.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.0278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[C]=O(1149)',
    structure = SMILES('[C]=O'),
    E0 = (440.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.66064,-0.00539267,9.3647e-06,-6.04676e-09,1.10218e-12,52863.3,2.60381], Tmin=(100,'K'), Tmax=(2084.48,'K')), NASAPolynomial(coeffs=[9.43361,-0.00191483,-2.23152e-06,5.70335e-10,-4.024e-14,48128.1,-30.5142], Tmin=(2084.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.031,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(CdCdJ2_triplet)"""),
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
    label = 'O=[C]C[C]=CO(28268)',
    structure = SMILES('O=[C]C[C]=CO'),
    E0 = (95.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1855,455,950,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.873885,'amu*angstrom^2'), symmetry=1, barrier=(20.0923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873934,'amu*angstrom^2'), symmetry=1, barrier=(20.0935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874179,'amu*angstrom^2'), symmetry=1, barrier=(20.0991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4357,0.0445082,-1.88996e-05,-1.7366e-08,1.21263e-11,11633.7,22.0205], Tmin=(100,'K'), Tmax=(972.355,'K')), NASAPolynomial(coeffs=[17.7368,0.00539319,-1.66531e-06,3.72394e-10,-3.30863e-14,7142.59,-62.9584], Tmin=(972.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=CO)C[C]=O(28269)',
    structure = SMILES('[CH]C(=CO)C[C]=O'),
    E0 = (189.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.98681,'amu*angstrom^2'), symmetry=1, barrier=(45.6806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98069,'amu*angstrom^2'), symmetry=1, barrier=(45.5399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98857,'amu*angstrom^2'), symmetry=1, barrier=(45.7212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98961,'amu*angstrom^2'), symmetry=1, barrier=(45.745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633672,0.0597233,-2.38227e-05,-2.08596e-08,1.43219e-11,22944.7,26.2912], Tmin=(100,'K'), Tmax=(978.353,'K')), NASAPolynomial(coeffs=[19.7054,0.0147276,-5.39921e-06,1.04109e-09,-7.87278e-14,17634.6,-73.3677], Tmin=(978.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(C=O)C[C]=O(14218)',
    structure = SMILES('C=C(C=O)C[C]=O'),
    E0 = (-63.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,462.726],'cm^-1')),
        HinderedRotor(inertia=(0.0939703,'amu*angstrom^2'), symmetry=1, barrier=(14.2806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939194,'amu*angstrom^2'), symmetry=1, barrier=(14.2797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939375,'amu*angstrom^2'), symmetry=1, barrier=(14.2795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80796,0.0476568,-3.10014e-05,8.43931e-09,-8.69774e-13,-7616.43,23.0065], Tmin=(100,'K'), Tmax=(2277.94,'K')), NASAPolynomial(coeffs=[22.7653,0.0108564,-6.76873e-06,1.34735e-09,-9.14449e-14,-17164.3,-95.3539], Tmin=(2277.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C(C[O])C[C]=O(9750)',
    structure = SMILES('C=C(C[O])C[C]=O'),
    E0 = (100.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(1.12365,'amu*angstrom^2'), symmetry=1, barrier=(25.8349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511152,'amu*angstrom^2'), symmetry=1, barrier=(25.861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1246,'amu*angstrom^2'), symmetry=1, barrier=(25.8568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48371,0.0446178,-2.50125e-05,5.8223e-09,-5.10003e-13,12110,23.2467], Tmin=(100,'K'), Tmax=(2625.4,'K')), NASAPolynomial(coeffs=[22.1236,0.0146948,-7.91624e-06,1.48105e-09,-9.66137e-14,1797.51,-90.4618], Tmin=(2625.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C[C]=O)CO(28270)',
    structure = SMILES('[CH2]C(=C[C]=O)CO'),
    E0 = (2.56546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.690654,'amu*angstrom^2'), symmetry=1, barrier=(15.8795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585584,'amu*angstrom^2'), symmetry=1, barrier=(15.8759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325014,'amu*angstrom^2'), symmetry=1, barrier=(7.47271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.641436,'amu*angstrom^2'), symmetry=1, barrier=(15.8706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830201,0.0757056,-0.000102552,7.87908e-08,-2.49449e-11,417.196,25.9493], Tmin=(100,'K'), Tmax=(765.529,'K')), NASAPolynomial(coeffs=[9.32136,0.0313368,-1.56124e-05,3.07669e-09,-2.182e-13,-882.816,-12.7467], Tmin=(765.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.56546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(CO)C[C]=O(28271)',
    structure = SMILES('[CH]=C(CO)C[C]=O'),
    E0 = (121.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3120,650,792.5,1650,350,440,435,1725,395.94],'cm^-1')),
        HinderedRotor(inertia=(0.00107219,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995496,'amu*angstrom^2'), symmetry=1, barrier=(11.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995584,'amu*angstrom^2'), symmetry=1, barrier=(11.0473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0993812,'amu*angstrom^2'), symmetry=1, barrier=(11.0518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53349,0.0564015,-4.70088e-05,1.98742e-08,-3.48631e-12,14727.1,26.2497], Tmin=(100,'K'), Tmax=(1309.56,'K')), NASAPolynomial(coeffs=[11.2478,0.0267298,-1.30225e-05,2.57267e-09,-1.83412e-13,12182.8,-23.236], Tmin=(1309.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C(=CO)CC=O(28272)',
    structure = SMILES('[CH]C(=CO)CC=O'),
    E0 = (29.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.00868,'amu*angstrom^2'), symmetry=1, barrier=(46.1834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00475,'amu*angstrom^2'), symmetry=1, barrier=(46.0932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00281,'amu*angstrom^2'), symmetry=1, barrier=(46.0486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00413,'amu*angstrom^2'), symmetry=1, barrier=(46.0789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688501,0.0569587,-9.41859e-06,-3.43937e-08,1.82833e-11,3705.18,25.2764], Tmin=(100,'K'), Tmax=(991.907,'K')), NASAPolynomial(coeffs=[19.0893,0.0192577,-7.60663e-06,1.48942e-09,-1.11607e-13,-1740.92,-72.3993], Tmin=(991.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C1C[C]([CH]O)C1(28273)',
    structure = SMILES('O=C1C[C]([CH]O)C1'),
    E0 = (24.2931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64517,0.043603,-1.30275e-05,-1.01622e-08,5.81168e-12,3013.49,23.6328], Tmin=(100,'K'), Tmax=(1083.35,'K')), NASAPolynomial(coeffs=[10.8111,0.0257851,-1.05451e-05,1.96422e-09,-1.37559e-13,87.1323,-25.6613], Tmin=(1083.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.2931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C=C(C=C=O)CO(28274)',
    structure = SMILES('C=C(C=C=O)CO'),
    E0 = (-166.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00028,0.070866,-9.84681e-05,7.99218e-08,-2.63934e-11,-19955.3,25.3702], Tmin=(100,'K'), Tmax=(808.238,'K')), NASAPolynomial(coeffs=[8.28089,0.0299996,-1.36525e-05,2.56198e-09,-1.75729e-13,-20974.3,-7.22769], Tmin=(808.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(C=O)CC=O(14231)',
    structure = SMILES('C=C(C=O)CC=O'),
    E0 = (-223.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68558,0.0467632,-2.22032e-05,9.6247e-10,9.90228e-13,-26847.7,22.6434], Tmin=(100,'K'), Tmax=(1552.48,'K')), NASAPolynomial(coeffs=[16.2483,0.0229833,-1.25038e-05,2.49863e-09,-1.75239e-13,-33025.3,-59.3518], Tmin=(1552.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=[C]C(O)C[C]=O(14295)',
    structure = SMILES('C=[C]C(O)C[C]=O'),
    E0 = (101.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1855,455,950,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,243.441,243.605],'cm^-1')),
        HinderedRotor(inertia=(0.00285574,'amu*angstrom^2'), symmetry=1, barrier=(0.119631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193957,'amu*angstrom^2'), symmetry=1, barrier=(8.10241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192572,'amu*angstrom^2'), symmetry=1, barrier=(8.09824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191549,'amu*angstrom^2'), symmetry=1, barrier=(8.09734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4010.29,'J/mol'), sigma=(6.50189,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=626.40 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17162,0.0659362,-7.8981e-05,5.36055e-08,-1.49744e-11,12353.6,29.0775], Tmin=(100,'K'), Tmax=(864.703,'K')), NASAPolynomial(coeffs=[9.44945,0.0276434,-1.25533e-05,2.39036e-09,-1.66946e-13,10922.1,-9.65477], Tmin=(864.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1CC(=O)C1O(28275)',
    structure = SMILES('C=C1CC(=O)C1O'),
    E0 = (-176.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02533,0.0221317,6.64578e-05,-1.01708e-07,3.92725e-11,-21153.1,21.6335], Tmin=(100,'K'), Tmax=(996.213,'K')), NASAPolynomial(coeffs=[16.2086,0.0183841,-8.00472e-06,1.72852e-09,-1.37553e-13,-26619,-59.9887], Tmin=(996.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=[C]C[C]=O(4742)',
    structure = SMILES('C=[C]C[C]=O'),
    E0 = (304.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1855,455,950,2950,3100,1380,975,1025,1650,471.004],'cm^-1')),
        HinderedRotor(inertia=(0.0843649,'amu*angstrom^2'), symmetry=1, barrier=(13.3205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842745,'amu*angstrom^2'), symmetry=1, barrier=(13.3224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49518,0.0282115,-1.11739e-05,-2.4489e-09,1.75756e-12,36700.7,18.3757], Tmin=(100,'K'), Tmax=(1305.79,'K')), NASAPolynomial(coeffs=[9.9248,0.0152979,-7.64936e-06,1.52571e-09,-1.08865e-13,33921,-22.6642], Tmin=(1305.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    E0 = (-29.5266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (169.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (158.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (172.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (166.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (82.3989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (173.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (225.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (94.5655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (129.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (65.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (426.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (323.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (347.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (331.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (422.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (201.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (96.7089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (33.8735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (346.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (315.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (216.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (235.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-21.2423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (361.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (114.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (535.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (511.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (401.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (168.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (212.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (167.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (266.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (74.0066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (96.7089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (33.8735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-4.55339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (196.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-21.2423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (544.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['C=C=O(598)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C1([CH]O)CC1=O(28256)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(199.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C(C=C=O)=CO(28257)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-OneDeH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C=O(598)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0561524,'m^3/(mol*s)'), n=2.47384, Ea=(45.2178,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C]=O(601)', 'C=C=CO(12571)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C(C=C[O])=CO(28258)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['CC([CH]O)=C[C]=O(28259)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['CC(=[C]O)C[C]=O(28260)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=[C]O)CC=O(28261)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['CC(=C[O])C[C]=O(14220)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(493835,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;XH_out] for rate rule [R4H_SDS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C(=C[O])CC=O(14224)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(D)(132)', '[CH]C(=C)C[C]=O(17871)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH2]C(=C[O])C[C]=O(14226)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(601)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2]C([CH]O)=C[C]=O(28262)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH2]C(=[C]O)C[C]=O(28263)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['O=[C]C[C]1CC1O(28264)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2][C]1CC(=O)C1O(28265)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['CC(C=C=O)=CO(28266)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C](C[O])C[C]=O(9762)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH]O)[CH][C]=O(14234)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC[C]=CO(14269)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['O=C1CC(=CO)C1(28267)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[C-]#[O+](374)', '[CH2]C([CH2])=CO(17684)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.52739e+06,'m^3/(mol*s)'), n=0.364815, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm
Ea raised from -181.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C(C=O)C[C]=O(12770)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[C]=O(1149)', '[CH2]C([CH2])=CO(17684)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(T)(28)', 'O=[C]C[C]=CO(28268)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(8)', '[CH]C(=CO)C[C]=O(28269)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(8)', 'C=C(C=O)C[C]=O(14218)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C(C[O])C[C]=O(9750)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=C[C]=O)CO(28270)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(400.242,'s^-1'), n=3.11383, Ea=(164.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;Cs_H_out_H/NonDeO] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_1H] for rate rule [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(CO)C[C]=O(28271)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(=CO)CC=O(28272)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['O=C1C[C]([CH]O)C1(28273)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['C=C(C=C=O)CO(28274)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['C=C(C=O)CC=O(14231)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(O)C[C]=O(14295)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['C=C1CC(=O)C1O(28275)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]C[C]=O(4742)', '[CH]O(5471)'],
    products = ['[CH2]C(=CO)C[C]=O(14236)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

network(
    label = '5114',
    isomers = [
        '[CH2]C(=CO)C[C]=O(14236)',
    ],
    reactants = [
        ('C=C=O(598)', 'C=C=CO(12571)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '5114',
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


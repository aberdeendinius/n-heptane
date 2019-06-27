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
    label = '[CH2]C1([O])CC1=CO(28311)',
    structure = SMILES('[CH2]C1([O])CC1=CO'),
    E0 = (199.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695759,0.0560246,-9.6232e-06,-4.46797e-08,2.62734e-11,24077.2,25.0873], Tmin=(100,'K'), Tmax=(917.677,'K')), NASAPolynomial(coeffs=[22.2422,0.00587837,7.97451e-07,-2.73511e-10,1.60879e-14,18279.6,-87.0526], Tmin=(917.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1([CH]O)OC1=C(28294)',
    structure = SMILES('[CH2]C1([CH]O)OC1=C'),
    E0 = (164.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.704672,0.085685,-0.000100354,5.4564e-08,-1.09053e-11,19924.8,25.9253], Tmin=(100,'K'), Tmax=(1422.84,'K')), NASAPolynomial(coeffs=[23.6454,0.00278445,2.27129e-06,-6.56152e-10,5.08744e-14,14457.8,-94.9986], Tmin=(1422.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCsJOH) + radical(CJC(C)OC)"""),
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
    label = 'C=C([O])C(C)=[C]O(28332)',
    structure = SMILES('C=C([O])C(C)=[C]O'),
    E0 = (11.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139851,0.078887,-9.05149e-05,5.01875e-08,-1.04241e-11,1505.86,27.1812], Tmin=(100,'K'), Tmax=(1329.73,'K')), NASAPolynomial(coeffs=[20.2342,0.00837587,-5.70601e-07,-1.22716e-10,1.52953e-14,-3097.1,-73.8516], Tmin=(1329.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([CH2])=CO(28333)',
    structure = SMILES('[CH]=C(O)C([CH2])=CO'),
    E0 = (32.2404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750],'cm^-1')),
        HinderedRotor(inertia=(1.16798,'amu*angstrom^2'), symmetry=1, barrier=(26.8543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16863,'amu*angstrom^2'), symmetry=1, barrier=(26.8692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16956,'amu*angstrom^2'), symmetry=1, barrier=(26.8905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16813,'amu*angstrom^2'), symmetry=1, barrier=(26.8577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96764,0.0957018,-0.000110386,5.66725e-08,-1.04017e-11,4123.18,29.8171], Tmin=(100,'K'), Tmax=(1627.91,'K')), NASAPolynomial(coeffs=[27.5646,-0.00405035,6.57971e-06,-1.48647e-09,1.05302e-13,-1889.51,-115.985], Tmin=(1627.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.2404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=[C]O)C(=C)O(28334)',
    structure = SMILES('[CH2]C(=[C]O)C(=C)O'),
    E0 = (24.8884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1685,370,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(1.03983,'amu*angstrom^2'), symmetry=1, barrier=(23.9076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03931,'amu*angstrom^2'), symmetry=1, barrier=(23.8957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03961,'amu*angstrom^2'), symmetry=1, barrier=(23.9027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04147,'amu*angstrom^2'), symmetry=1, barrier=(23.9454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15441,0.087855,-9.99809e-05,5.20255e-08,-9.82832e-12,3200.92,29.6121], Tmin=(100,'K'), Tmax=(1545.63,'K')), NASAPolynomial(coeffs=[24.977,0.000585884,3.77445e-06,-9.48672e-10,6.99519e-14,-2530.71,-100.245], Tmin=(1545.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.8884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C(C)=CO(28335)',
    structure = SMILES('[CH]=C([O])C(C)=CO'),
    E0 = (18.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(0.954039,'amu*angstrom^2'), symmetry=1, barrier=(21.9352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956844,'amu*angstrom^2'), symmetry=1, barrier=(21.9997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.953101,'amu*angstrom^2'), symmetry=1, barrier=(21.9137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.875888,0.0859393,-9.8619e-05,5.24294e-08,-1.01843e-11,2424.45,27.1007], Tmin=(100,'K'), Tmax=(1487.32,'K')), NASAPolynomial(coeffs=[23.8325,0.00239514,2.87655e-06,-7.91479e-10,6.03244e-14,-3034.78,-95.5559], Tmin=(1487.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C(C)=C[O](14400)',
    structure = SMILES('C=C([O])C(C)=C[O]'),
    E0 = (-87.0875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.38725,0.0751821,-7.66173e-05,3.73455e-08,-6.75342e-12,-10297.9,26.3693], Tmin=(100,'K'), Tmax=(1580.12,'K')), NASAPolynomial(coeffs=[21.0081,0.00785301,-2.02134e-07,-1.68594e-10,1.63062e-14,-15415.5,-81.4377], Tmin=(1580.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.0875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=C[O])C(=C)O(14402)',
    structure = SMILES('[CH2]C(=C[O])C(=C)O'),
    E0 = (-73.3931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24029,'amu*angstrom^2'), symmetry=1, barrier=(28.5167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2366,'amu*angstrom^2'), symmetry=1, barrier=(28.432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23687,'amu*angstrom^2'), symmetry=1, barrier=(28.4381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.172681,0.0658946,-2.34589e-05,-4.07186e-08,2.75213e-11,-8672.22,23.1281], Tmin=(100,'K'), Tmax=(907.998,'K')), NASAPolynomial(coeffs=[26.3498,-3.4751e-05,3.86626e-06,-8.77071e-10,5.82081e-14,-15461.9,-111.846], Tmin=(907.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.3931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=COJ)"""),
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
    label = '[CH]C(=C)C(=C)[O](17908)',
    structure = SMILES('[CH]C(=C)C(=C)[O]'),
    E0 = (350.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15646,'amu*angstrom^2'), symmetry=1, barrier=(49.5813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16216,'amu*angstrom^2'), symmetry=1, barrier=(49.7124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.45,'J/mol'), sigma=(6.08582,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.10 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15771,0.0556689,-3.86778e-05,5.32957e-09,3.75398e-12,42315.3,20.6503], Tmin=(100,'K'), Tmax=(941.054,'K')), NASAPolynomial(coeffs=[13.7286,0.0187804,-6.25033e-06,1.03926e-09,-6.93387e-14,39216.8,-43.1253], Tmin=(941.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(=C[O])C(=C)[O](14406)',
    structure = SMILES('[CH2]C(=C[O])C(=C)[O]'),
    E0 = (64.4117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30251,'amu*angstrom^2'), symmetry=1, barrier=(29.9473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2993,'amu*angstrom^2'), symmetry=1, barrier=(29.8735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.705456,0.0763061,-8.01805e-05,3.90496e-08,-6.92712e-12,7939.82,27.8011], Tmin=(100,'K'), Tmax=(1656.03,'K')), NASAPolynomial(coeffs=[21.9753,0.00356552,1.97173e-06,-5.7025e-10,4.25056e-14,2890.15,-85.6267], Tmin=(1656.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=[C]O)C(=C)[O](28336)',
    structure = SMILES('[CH2]C(=[C]O)C(=C)[O]'),
    E0 = (162.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.960158,'amu*angstrom^2'), symmetry=1, barrier=(22.0759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960296,'amu*angstrom^2'), symmetry=1, barrier=(22.0791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960498,'amu*angstrom^2'), symmetry=1, barrier=(22.0837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.388544,0.079299,-9.20287e-05,4.97626e-08,-9.88196e-12,19740.3,28.3555], Tmin=(100,'K'), Tmax=(1439.85,'K')), NASAPolynomial(coeffs=[21.9626,0.00308601,2.0778e-06,-6.20513e-10,4.85587e-14,14767.5,-82.5412], Tmin=(1439.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C([CH2])=CO(28337)',
    structure = SMILES('[CH]=C([O])C([CH2])=CO'),
    E0 = (170.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(1.15249,'amu*angstrom^2'), symmetry=1, barrier=(26.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14833,'amu*angstrom^2'), symmetry=1, barrier=(26.4023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14962,'amu*angstrom^2'), symmetry=1, barrier=(26.4321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16454,0.0867695,-0.000101372,5.3335e-08,-1.01046e-11,20660.8,28.4223], Tmin=(100,'K'), Tmax=(1557.44,'K')), NASAPolynomial(coeffs=[24.9939,-0.00212208,5.14859e-06,-1.2112e-09,8.77374e-14,15145.5,-100.914], Tmin=(1557.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[C]1CC1O(28338)',
    structure = SMILES('C=C([O])[C]1CC1O'),
    E0 = (17.6722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51478,0.0372001,2.79917e-05,-7.20213e-08,3.29304e-11,2231.02,24.6624], Tmin=(100,'K'), Tmax=(934.593,'K')), NASAPolynomial(coeffs=[17.8865,0.0115529,-2.14271e-06,3.3252e-10,-2.83849e-14,-2769.25,-63.5937], Tmin=(934.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.6722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJ(C)CO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([CH]O)=C1CO1(28339)',
    structure = SMILES('[CH2]C([CH]O)=C1CO1'),
    E0 = (64.8629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06601,0.043616,2.59451e-05,-8.05753e-08,3.84431e-11,7926.17,23.199], Tmin=(100,'K'), Tmax=(928.274,'K')), NASAPolynomial(coeffs=[22.4196,0.00534788,9.34338e-07,-2.40334e-10,9.50113e-15,1646.15,-90.7032], Tmin=(928.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.8629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[O]C1CCC=1[CH]O(28340)',
    structure = SMILES('[O]C1CCC=1[CH]O'),
    E0 = (-7.1819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50125,0.0357109,3.64201e-05,-8.26744e-08,3.70175e-11,-756.124,24.9361], Tmin=(100,'K'), Tmax=(936.332,'K')), NASAPolynomial(coeffs=[18.8291,0.0106012,-1.71564e-06,2.71275e-10,-2.56979e-14,-6145.28,-68.9711], Tmin=(936.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.1819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=CCJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1OC(O)C=1[CH2](28341)',
    structure = SMILES('[CH2]C1OC(O)C=1[CH2]'),
    E0 = (27.4661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.805467,0.0532571,-4.44114e-06,-4.63034e-08,2.5399e-11,3434.05,22.6571], Tmin=(100,'K'), Tmax=(937.514,'K')), NASAPolynomial(coeffs=[21.4244,0.00822035,-1.07994e-06,1.56582e-10,-1.66826e-14,-2318.94,-85.55], Tmin=(937.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.4661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])=C([CH2])C[O](9827)',
    structure = SMILES('[CH2]C([O])=C([CH2])C[O]'),
    E0 = (237.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1016.49,1016.99],'cm^-1')),
        HinderedRotor(inertia=(0.112021,'amu*angstrom^2'), symmetry=1, barrier=(2.57558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111837,'amu*angstrom^2'), symmetry=1, barrier=(2.57136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111721,'amu*angstrom^2'), symmetry=1, barrier=(2.56868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08165,0.0669476,-7.93282e-05,5.43125e-08,-1.52735e-11,28640.1,27.3547], Tmin=(100,'K'), Tmax=(861.487,'K')), NASAPolynomial(coeffs=[9.41492,0.0282583,-1.1969e-05,2.19052e-09,-1.49181e-13,27204.2,-11.6068], Tmin=(861.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(O)CJ) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=[C]O)C([CH2])[O](28342)',
    structure = SMILES('[CH2]C(=[C]O)C([CH2])[O]'),
    E0 = (396.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,320.036,320.295],'cm^-1')),
        HinderedRotor(inertia=(0.293294,'amu*angstrom^2'), symmetry=1, barrier=(21.3431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293248,'amu*angstrom^2'), symmetry=1, barrier=(21.3419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293326,'amu*angstrom^2'), symmetry=1, barrier=(21.3408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164424,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313512,0.0781262,-8.35159e-05,4.24014e-08,-8.05547e-12,47799.5,33.3], Tmin=(100,'K'), Tmax=(1455.2,'K')), NASAPolynomial(coeffs=[21.9965,0.00729085,-6.96541e-07,-3.15596e-11,5.98266e-15,42313.3,-79.2425], Tmin=(1455.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CJCO) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])[CH]O(14416)',
    structure = SMILES('[CH]=C([O])C([CH2])[CH]O'),
    E0 = (347.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,323.459,323.46],'cm^-1')),
        HinderedRotor(inertia=(0.00161112,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00161115,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128675,'amu*angstrom^2'), symmetry=1, barrier=(9.55395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128678,'amu*angstrom^2'), symmetry=1, barrier=(9.55384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482846,0.0786773,-9.88423e-05,5.88739e-08,-1.14794e-11,41934.9,30.7124], Tmin=(100,'K'), Tmax=(768.502,'K')), NASAPolynomial(coeffs=[14.8639,0.0173504,-5.54127e-06,8.37961e-10,-4.99517e-14,39325.1,-37.4799], Tmin=(768.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl) + radical(CCsJOH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([C]O)=C(C)[O](28343)',
    structure = SMILES('[CH2]C([C]O)=C(C)[O]'),
    E0 = (313.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,345.631,345.875,345.971],'cm^-1')),
        HinderedRotor(inertia=(0.00140767,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136301,'amu*angstrom^2'), symmetry=1, barrier=(11.5801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136443,'amu*angstrom^2'), symmetry=1, barrier=(11.5827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135976,'amu*angstrom^2'), symmetry=1, barrier=(11.5744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534565,0.0752663,-8.99139e-05,5.49947e-08,-1.31916e-11,37848,26.4952], Tmin=(100,'K'), Tmax=(1023.72,'K')), NASAPolynomial(coeffs=[14.9217,0.0190508,-7.54397e-06,1.3534e-09,-9.195e-14,34902.3,-43.2516], Tmin=(1023.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C([O])=C([CH2])CO(28344)',
    structure = SMILES('[CH]C([O])=C([CH2])CO'),
    E0 = (223.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675629,0.0701193,-6.75415e-05,3.53631e-08,-7.50518e-12,26983.5,28.3632], Tmin=(100,'K'), Tmax=(1135.48,'K')), NASAPolynomial(coeffs=[12.9483,0.0268859,-1.04291e-05,1.83123e-09,-1.22459e-13,24196.4,-32.405], Tmin=(1135.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(AllylJ2_triplet) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=C[O])C([CH2])[O](14412)',
    structure = SMILES('[CH2]C(=C[O])C([CH2])[O]'),
    E0 = (297.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,444.615,444.732,445.031],'cm^-1')),
        HinderedRotor(inertia=(0.165217,'amu*angstrom^2'), symmetry=1, barrier=(23.2162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165524,'amu*angstrom^2'), symmetry=1, barrier=(23.2102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165305,'amu*angstrom^2'), symmetry=1, barrier=(23.2126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509394,0.0621704,-2.83882e-05,-2.18082e-08,1.66898e-11,35948.4,28.6222], Tmin=(100,'K'), Tmax=(942.148,'K')), NASAPolynomial(coeffs=[21.53,0.00921499,-1.8551e-06,3.00587e-10,-2.5321e-14,30376.8,-80.0851], Tmin=(942.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CJCO) + radical(CC(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH][O])=C(C)[O](28345)',
    structure = SMILES('[CH2]C([CH][O])=C(C)[O]'),
    E0 = (195.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,389.499,390.103,390.658],'cm^-1')),
        HinderedRotor(inertia=(0.00110524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111033,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00110224,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00491,0.0591228,-4.9203e-05,2.11003e-08,-3.63192e-12,23645,27.5365], Tmin=(100,'K'), Tmax=(1387.96,'K')), NASAPolynomial(coeffs=[14.1084,0.0213594,-8.391e-06,1.49735e-09,-1.01019e-13,20007.6,-39.976], Tmin=(1387.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.56428,0.0758759,-7.90326e-05,3.87517e-08,-6.96524e-12,5122.04,31.1132], Tmin=(100,'K'), Tmax=(1620.33,'K')), NASAPolynomial(coeffs=[21.4332,0.00536805,1.23983e-06,-4.47392e-10,3.50195e-14,120.562,-79.0647], Tmin=(1620.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.0451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OCC1=CO(28283)',
    structure = SMILES('C=C1OCC1=CO'),
    E0 = (-165.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.880308,0.0456915,2.88486e-05,-9.02484e-08,4.37695e-11,-19809.6,18.361], Tmin=(100,'K'), Tmax=(917.983,'K')), NASAPolynomial(coeffs=[24.8046,0.00132924,3.48413e-06,-7.64082e-10,4.64529e-14,-26725.2,-108.756], Tmin=(917.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = '[CH]C(=C)C(=C)OO(20661)',
    structure = SMILES('[CH]C(=C)C(=C)OO'),
    E0 = (337.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,217.836,217.836,217.836,217.836],'cm^-1')),
        HinderedRotor(inertia=(1.48426,'amu*angstrom^2'), symmetry=1, barrier=(49.9796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48425,'amu*angstrom^2'), symmetry=1, barrier=(49.9796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48424,'amu*angstrom^2'), symmetry=1, barrier=(49.9796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48425,'amu*angstrom^2'), symmetry=1, barrier=(49.9796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725505,0.0722406,-7.04379e-05,3.72786e-08,-8.11007e-12,40723.7,24.7482], Tmin=(100,'K'), Tmax=(1096.53,'K')), NASAPolynomial(coeffs=[12.2055,0.0303635,-1.31525e-05,2.45054e-09,-1.69629e-13,38206.1,-31.694], Tmin=(1096.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    label = 'C=C([O])[C]=CO(28346)',
    structure = SMILES('C=C([O])[C]=CO'),
    E0 = (8.03594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29004,'amu*angstrom^2'), symmetry=1, barrier=(29.6606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28562,'amu*angstrom^2'), symmetry=1, barrier=(29.5589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.102175,0.0698031,-8.38058e-05,4.56626e-08,-8.91755e-12,1131.82,22.9119], Tmin=(100,'K'), Tmax=(1525.67,'K')), NASAPolynomial(coeffs=[19.8181,-0.00133889,4.73582e-06,-1.15325e-09,8.5368e-14,-2745.13,-74.3923], Tmin=(1525.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.03594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C([C]=C)=CO(28222)',
    structure = SMILES('[CH2]C([C]=C)=CO'),
    E0 = (198.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.44396,'amu*angstrom^2'), symmetry=1, barrier=(33.1994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44198,'amu*angstrom^2'), symmetry=1, barrier=(33.154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.443,'amu*angstrom^2'), symmetry=1, barrier=(33.1775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.780012,0.0740694,-7.72484e-05,3.75677e-08,-6.5926e-12,24016.9,25.5169], Tmin=(100,'K'), Tmax=(1711.77,'K')), NASAPolynomial(coeffs=[20.126,0.00408743,2.59141e-06,-7.3779e-10,5.49892e-14,19955.3,-77.5378], Tmin=(1711.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C(=CO)C(=C)[O](28347)',
    structure = SMILES('[CH]C(=CO)C(=C)[O]'),
    E0 = (142.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.98945,'amu*angstrom^2'), symmetry=1, barrier=(45.7415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98966,'amu*angstrom^2'), symmetry=1, barrier=(45.7461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98978,'amu*angstrom^2'), symmetry=1, barrier=(45.749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06586,0.0852754,-9.1073e-05,4.57823e-08,-8.42975e-12,17299.7,28.5001], Tmin=(100,'K'), Tmax=(1585.4,'K')), NASAPolynomial(coeffs=[23.0281,0.0059078,1.59658e-06,-5.76839e-10,4.60921e-14,11994.8,-91.4795], Tmin=(1585.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1([O])C(=C)C1O(28320)',
    structure = SMILES('[CH2]C1([O])C(=C)C1O'),
    E0 = (235.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02709,0.0549503,-2.84615e-05,-7.97311e-09,8.45389e-12,28387.2,25.9186], Tmin=(100,'K'), Tmax=(980.108,'K')), NASAPolynomial(coeffs=[16.1256,0.0167676,-5.89423e-06,1.07482e-09,-7.73369e-14,24301.8,-52.3623], Tmin=(980.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C(=C)C=O(14399)',
    structure = SMILES('C=C([O])C(=C)C=O'),
    E0 = (-96.1204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.921994,'amu*angstrom^2'), symmetry=1, barrier=(21.1985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920976,'amu*angstrom^2'), symmetry=1, barrier=(21.175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510012,0.067351,-6.87242e-05,3.39519e-08,-6.46547e-12,-11427,22.5093], Tmin=(100,'K'), Tmax=(1293.46,'K')), NASAPolynomial(coeffs=[18.4591,0.0118447,-4.35555e-06,7.75925e-10,-5.32954e-14,-16070.3,-68.7034], Tmin=(1293.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.1204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(=C)C[O](9813)',
    structure = SMILES('C=C([O])C(=C)C[O]'),
    E0 = (52.8118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,299.045,299.265,299.335,2642.21],'cm^-1')),
        HinderedRotor(inertia=(0.22376,'amu*angstrom^2'), symmetry=1, barrier=(14.238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224151,'amu*angstrom^2'), symmetry=1, barrier=(14.2361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13901,0.0654443,-7.4141e-05,4.7983e-08,-1.2779e-11,6452.63,24.2329], Tmin=(100,'K'), Tmax=(906.369,'K')), NASAPolynomial(coeffs=[9.66234,0.0278303,-1.18937e-05,2.19951e-09,-1.51165e-13,4907.52,-16.0496], Tmin=(906.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.8118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(CO)C(=C)[O](28348)',
    structure = SMILES('[CH]=C(CO)C(=C)[O]'),
    E0 = (74.2029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,251.093,251.129],'cm^-1')),
        HinderedRotor(inertia=(0.255418,'amu*angstrom^2'), symmetry=1, barrier=(11.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255423,'amu*angstrom^2'), symmetry=1, barrier=(11.4286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255472,'amu*angstrom^2'), symmetry=1, barrier=(11.4289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498762,0.074265,-8.81002e-05,5.36979e-08,-1.27538e-11,9053.01,26.0637], Tmin=(100,'K'), Tmax=(1038.43,'K')), NASAPolynomial(coeffs=[15.2545,0.0174258,-5.99575e-06,9.86617e-10,-6.35159e-14,5988.49,-45.6803], Tmin=(1038.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.2029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=CO)C(=C)O(28349)',
    structure = SMILES('[CH]C(=CO)C(=C)O'),
    E0 = (4.3295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.89146,'amu*angstrom^2'), symmetry=1, barrier=(43.4884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89096,'amu*angstrom^2'), symmetry=1, barrier=(43.4768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89049,'amu*angstrom^2'), symmetry=1, barrier=(43.4662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88941,'amu*angstrom^2'), symmetry=1, barrier=(43.4413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88266,0.0943418,-0.000100447,4.94648e-08,-8.83287e-12,762.794,29.9462], Tmin=(100,'K'), Tmax=(1660.43,'K')), NASAPolynomial(coeffs=[25.6023,0.0039539,3.04835e-06,-8.57758e-10,6.41628e-14,-5031.84,-106.554], Tmin=(1660.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.3295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C(=C)CO(28350)',
    structure = SMILES('[CH]=C([O])C(=C)CO'),
    E0 = (74.2029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,251.093,251.129],'cm^-1')),
        HinderedRotor(inertia=(0.255418,'amu*angstrom^2'), symmetry=1, barrier=(11.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255423,'amu*angstrom^2'), symmetry=1, barrier=(11.4286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255472,'amu*angstrom^2'), symmetry=1, barrier=(11.4289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498762,0.074265,-8.81002e-05,5.36979e-08,-1.27538e-11,9053.01,26.0637], Tmin=(100,'K'), Tmax=(1038.43,'K')), NASAPolynomial(coeffs=[15.2545,0.0174258,-5.99575e-06,9.86617e-10,-6.35159e-14,5988.49,-45.6803], Tmin=(1038.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.2029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1=C([O])CC1O(28351)',
    structure = SMILES('[CH2]C1=C([O])CC1O'),
    E0 = (7.35123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23389,0.0432421,1.6521e-05,-6.35282e-08,3.0579e-11,999.914,24.3982], Tmin=(100,'K'), Tmax=(936.559,'K')), NASAPolynomial(coeffs=[19.4746,0.00995283,-1.62002e-06,2.50191e-10,-2.32802e-14,-4373.52,-72.8536], Tmin=(936.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.35123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1OCC=1[CH]O(28352)',
    structure = SMILES('[CH2]C1OCC=1[CH]O'),
    E0 = (47.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19281,0.0405714,3.30212e-05,-8.63734e-08,4.01095e-11,5779.31,23.9241], Tmin=(100,'K'), Tmax=(929.364,'K')), NASAPolynomial(coeffs=[21.7411,0.00643482,4.71196e-07,-1.51934e-10,3.19676e-15,-385.22,-86.3215], Tmin=(929.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(O)C(=C)C=O(14411)',
    structure = SMILES('C=C(O)C(=C)C=O'),
    E0 = (-233.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325514,0.0692779,-5.45287e-05,8.93633e-09,4.60925e-12,-27992.1,21.6647], Tmin=(100,'K'), Tmax=(977.967,'K')), NASAPolynomial(coeffs=[20.849,0.0106793,-3.52389e-06,6.66569e-10,-5.08648e-14,-33218.3,-83.0885], Tmin=(977.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]C(=CO)C([CH2])[O](23305)',
    structure = SMILES('[CH]C(=CO)C([CH2])[O]'),
    E0 = (375.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.104473,0.0714777,-3.96389e-05,-1.5715e-08,1.59161e-11,45310.6,29.4941], Tmin=(100,'K'), Tmax=(924.894,'K')), NASAPolynomial(coeffs=[22.4254,0.0116865,-2.25916e-06,2.94231e-10,-2.13986e-14,39610.1,-84.9448], Tmin=(924.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CJCO) + radical(AllylJ2_triplet) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]C([CH]O)=C(C)[O](23303)',
    structure = SMILES('[CH]C([CH]O)=C(C)[O]'),
    E0 = (189.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734299,0.0633331,-4.24096e-05,6.49637e-09,3.12392e-12,22872.5,28.4686], Tmin=(100,'K'), Tmax=(989.998,'K')), NASAPolynomial(coeffs=[14.8665,0.02364,-8.6425e-06,1.51787e-09,-1.03824e-13,19221.4,-43.8773], Tmin=(989.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=CCJO) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([O])=C(C)[CH]O(28353)',
    structure = SMILES('[CH]C([O])=C(C)[CH]O'),
    E0 = (189.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734299,0.0633331,-4.24096e-05,6.49637e-09,3.12392e-12,22872.5,28.4686], Tmin=(100,'K'), Tmax=(989.998,'K')), NASAPolynomial(coeffs=[14.8665,0.02364,-8.6425e-06,1.51787e-09,-1.03824e-13,19221.4,-43.8773], Tmin=(989.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=CCJO) + radical(AllylJ2_triplet) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C(O)C(=C)[O](14497)',
    structure = SMILES('C=[C]C(O)C(=C)[O]'),
    E0 = (72.0046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,350,440,435,1725,1380,1390,370,380,2900,435,418.502,418.788,419.084],'cm^-1')),
        HinderedRotor(inertia=(0.0748184,'amu*angstrom^2'), symmetry=1, barrier=(9.33062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751913,'amu*angstrom^2'), symmetry=1, barrier=(9.33013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749943,'amu*angstrom^2'), symmetry=1, barrier=(9.32591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4228.96,'J/mol'), sigma=(6.78198,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=660.55 K, Pc=30.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910647,0.0639176,-6.36516e-05,3.26891e-08,-6.65416e-12,8774.8,29.0222], Tmin=(100,'K'), Tmax=(1193.58,'K')), NASAPolynomial(coeffs=[14.2257,0.0192951,-7.57309e-06,1.36659e-09,-9.3486e-14,5596.3,-37.5713], Tmin=(1193.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.0046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC(O)C1=C(28293)',
    structure = SMILES('C=C1OC(O)C1=C'),
    E0 = (-163.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10621,0.049597,-6.93062e-06,-3.34333e-08,1.8026e-11,-19601.5,18.4114], Tmin=(100,'K'), Tmax=(968.551,'K')), NASAPolynomial(coeffs=[17.6605,0.0147637,-4.91867e-06,9.29083e-10,-7.05339e-14,-24381.1,-69.0445], Tmin=(968.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'C=[C]C(=C)[O](4767)',
    structure = SMILES('C=[C]C(=C)[O]'),
    E0 = (216.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44423,'amu*angstrom^2'), symmetry=1, barrier=(33.2056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65189,0.0455512,-4.93436e-05,2.74298e-08,-5.80388e-12,26168.2,16.7593], Tmin=(100,'K'), Tmax=(1302.66,'K')), NASAPolynomial(coeffs=[11.8538,0.00923597,-1.78226e-06,1.49149e-10,-4.09118e-15,23933.6,-33.5317], Tmin=(1302.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
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
    label = '[CH2]C(=[C]O)C(C)=O(28354)',
    structure = SMILES('[CH2]C(=[C]O)C(C)=O'),
    E0 = (17.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.645096,'amu*angstrom^2'), symmetry=1, barrier=(14.832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00565799,'amu*angstrom^2'), symmetry=1, barrier=(2.35082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0354814,'amu*angstrom^2'), symmetry=1, barrier=(14.818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646123,'amu*angstrom^2'), symmetry=1, barrier=(14.8556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88275,0.0627597,-5.86195e-05,2.78787e-08,-5.25512e-12,2217.74,26.4426], Tmin=(100,'K'), Tmax=(1283.94,'K')), NASAPolynomial(coeffs=[14.9089,0.0190628,-7.56959e-06,1.37185e-09,-9.39221e-14,-1384.01,-44.7309], Tmin=(1283.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2]C(=C[O])C(C)=O(14420)',
    structure = SMILES('[CH2]C(=C[O])C(C)=O'),
    E0 = (-80.8156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08232,0.0540748,-2.86103e-05,-4.0819e-09,5.83265e-12,-9606.05,24.0069], Tmin=(100,'K'), Tmax=(1032.9,'K')), NASAPolynomial(coeffs=[15.2211,0.0195553,-7.8652e-06,1.49451e-09,-1.07535e-13,-13606.2,-49.8873], Tmin=(1032.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.8156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + radical(C=C(C=O)CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)=C([CH2])[C]O(28355)',
    structure = SMILES('[CH2]C(O)=C([CH2])[C]O'),
    E0 = (334.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,305.502,305.915],'cm^-1')),
        HinderedRotor(inertia=(0.26382,'amu*angstrom^2'), symmetry=1, barrier=(17.4906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00180337,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263599,'amu*angstrom^2'), symmetry=1, barrier=(17.4894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264904,'amu*angstrom^2'), symmetry=1, barrier=(17.4893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263776,'amu*angstrom^2'), symmetry=1, barrier=(17.4888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261173,0.0845016,-0.000102441,5.93503e-08,-1.2984e-11,40422.9,28.5335], Tmin=(100,'K'), Tmax=(1207.51,'K')), NASAPolynomial(coeffs=[21.2858,0.00810973,-1.31517e-06,7.88043e-11,-4.30372e-16,35584.9,-77.967], Tmin=(1207.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CH2_triplet) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(O)=C([CH2])[CH][O](28356)',
    structure = SMILES('[CH2]C(O)=C([CH2])[CH][O]'),
    E0 = (216.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,607.096,607.096],'cm^-1')),
        HinderedRotor(inertia=(0.0107162,'amu*angstrom^2'), symmetry=1, barrier=(2.80273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930585,'amu*angstrom^2'), symmetry=1, barrier=(24.3387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930583,'amu*angstrom^2'), symmetry=1, barrier=(24.3387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306133,'amu*angstrom^2'), symmetry=1, barrier=(80.0666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.682559,0.062799,-4.24856e-05,9.27524e-10,6.80726e-12,26199.5,27.8761], Tmin=(100,'K'), Tmax=(954.309,'K')), NASAPolynomial(coeffs=[17.8831,0.0146911,-4.57418e-06,7.837e-10,-5.54866e-14,21824.2,-60.0254], Tmin=(954.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=CCJO) + radical(CCOJ) + radical(Allyl_P)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4086.37,'J/mol'), sigma=(6.56031,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=638.28 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636097,0.057278,-1.55871e-05,-3.36923e-08,2.00748e-11,-3415,25.454], Tmin=(100,'K'), Tmax=(963.153,'K')), NASAPolynomial(coeffs=[22.0242,0.00863463,-2.40919e-06,5.0148e-10,-4.37789e-14,-9398.75,-86.6034], Tmin=(963.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Allyl_P)"""),
)

species(
    label = 'O=C1CCC1=CO(28245)',
    structure = SMILES('O=C1CCC1=CO'),
    E0 = (-220.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56023,0.0364248,2.8489e-05,-6.77184e-08,2.94768e-11,-26471.3,21.3945], Tmin=(100,'K'), Tmax=(965.909,'K')), NASAPolynomial(coeffs=[16.7964,0.0157387,-5.2462e-06,1.02146e-09,-7.974e-14,-31393,-61.8238], Tmin=(965.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C([C]=O)=CO(28357)',
    structure = SMILES('[CH2]C([C]=O)=CO'),
    E0 = (-13.5195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.05839,'amu*angstrom^2'), symmetry=1, barrier=(24.3344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05787,'amu*angstrom^2'), symmetry=1, barrier=(24.3225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05931,'amu*angstrom^2'), symmetry=1, barrier=(24.3555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06851,0.0541419,-4.0665e-05,1.0596e-09,6.80454e-12,-1511.07,20.2326], Tmin=(100,'K'), Tmax=(955.011,'K')), NASAPolynomial(coeffs=[18.8869,0.00410755,-7.11277e-07,1.37903e-10,-1.40243e-14,-6036.1,-70.7833], Tmin=(955.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.5195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH]C(=CO)C(C)=O(28358)',
    structure = SMILES('[CH]C(=CO)C(C)=O'),
    E0 = (-7.94525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,375,552.5,462.5,1710,350,440,435,1725,327.746,327.746,327.746,327.746],'cm^-1')),
        HinderedRotor(inertia=(0.649937,'amu*angstrom^2'), symmetry=1, barrier=(49.542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649943,'amu*angstrom^2'), symmetry=1, barrier=(49.542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649937,'amu*angstrom^2'), symmetry=1, barrier=(49.542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64994,'amu*angstrom^2'), symmetry=1, barrier=(49.542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708837,0.062918,-4.07409e-05,5.66657e-09,2.84693e-12,-828.897,25.5455], Tmin=(100,'K'), Tmax=(1035.85,'K')), NASAPolynomial(coeffs=[15.3527,0.0238035,-9.34502e-06,1.70833e-09,-1.19195e-13,-4797.94,-50.1329], Tmin=(1035.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.94525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(C=O)C(C)=O(14422)',
    structure = SMILES('C=C(C=O)C(C)=O'),
    E0 = (-253.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56701,0.0454831,-2.62535e-05,6.22057e-09,-5.54825e-13,-30470.1,18.7814], Tmin=(100,'K'), Tmax=(2538.2,'K')), NASAPolynomial(coeffs=[20.5932,0.0170753,-9.46533e-06,1.81111e-09,-1.20516e-13,-39621,-84.9754], Tmin=(2538.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[CH]C([CH]O)=C([CH2])O(23433)',
    structure = SMILES('[CH]C([CH]O)=C([CH2])O'),
    E0 = (210.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367352,0.0674119,-3.64474e-05,-1.37694e-08,1.41204e-11,25429,28.9755], Tmin=(100,'K'), Tmax=(923.401,'K')), NASAPolynomial(coeffs=[19.983,0.01483,-3.64597e-06,5.34843e-10,-3.65138e-14,20425.5,-71.5723], Tmin=(923.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=CCJO) + radical(AllylJ2_triplet) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1C(=O)CC1O(28255)',
    structure = SMILES('C=C1C(=O)CC1O'),
    E0 = (-184.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90253,0.0350968,1.11949e-05,-3.41099e-08,1.35279e-11,-22161.7,22.194], Tmin=(100,'K'), Tmax=(1070.78,'K')), NASAPolynomial(coeffs=[11.2595,0.0256909,-1.14171e-05,2.25015e-09,-1.63449e-13,-25630.1,-30.4266], Tmin=(1070.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    E0 = (-77.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (199.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (164.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (166.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (155.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (178.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (224.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (69.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (62.8546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (81.9849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (14.8893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (379.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (276.216,'kJ/mol'),
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
    E0 = (374.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (382.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (27.1265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (110.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (53.3938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (49.1846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (265.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (459.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (411.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (322.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (248.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (322.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (220.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (211.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-68.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (450.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (66.8158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (423.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (604.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (353.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (235.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (136.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (164.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (219.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (48.6381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (118.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (48.2686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (49.1846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-52.0777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (438.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (214.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (214.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (166.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-68.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (457.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (168.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (61.7745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (12.1728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (343.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (241.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (216.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (-68.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (402.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (36.3633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (-52.0777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (235.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (-68.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C=O(598)', 'C=C=CO(12571)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1([O])CC1=CO(28311)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(276.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 274.2 to 276.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1([CH]O)OC1=C(28294)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(241.974,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C]=O(601)', 'C=C=CO(12571)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=O(598)', '[CH2][C]=CO(18753)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(238.054,'m^3/(mol*s)'), n=1.55554, Ea=(28.5522,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C([O])C(C)=[C]O(28332)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(O)C([CH2])=CO(28333)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=[C]O)C(=C)O(28334)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C([O])C(C)=CO(28335)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C([O])C(C)=C[O](14400)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(493835,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;XH_out] for rate rule [R4H_SDS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=C[O])C(=C)O(14402)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(D)(132)', '[CH]C(=C)C(=C)[O](17908)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.05166e+07,'m^3/(mol*s)'), n=0.045, Ea=(0.1046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', '[CH2]C(=C[O])C(=C)[O](14406)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
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
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2]C(=[C]O)C(=C)[O](28336)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]=C([O])C([CH2])=CO(28337)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C([O])[C]1CC1O(28338)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(104.177,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_HNd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C([CH]O)=C1CO1(28339)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[O]C1CCC=1[CH]O(28340)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1OC(O)C=1[CH2](28341)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])=C([CH2])C[O](9827)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=[C]O)C([CH2])[O](28342)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C([O])C([CH2])[CH]O(14416)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([C]O)=C(C)[O](28343)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C([O])=C([CH2])CO(28344)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(=C[O])C([CH2])[O](14412)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH][O])=C(C)[O](28345)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C[C]=CO(14467)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C1OCC1=CO(28283)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=C)C(=C)OO(20661)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.95074e+10,'s^-1'), n=0, Ea=(112.549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OOH;Y_rad_out] for rate rule [R3OOH_DS;Cd_rad_out_H]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C(C=O)C(=C)[O](12762)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(28)', 'C=C([O])[C]=CO(28346)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(T)(63)', '[CH2]C([C]=C)=CO(28222)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(8)', '[CH]C(=CO)C(=C)[O](28347)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1([O])C(=C)C1O(28320)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.15385e+14,'s^-1'), n=-0.537569, Ea=(312.108,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 309.6 to 312.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', 'C=C([O])C(=C)C=O(14399)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2834 used for Od_CO-CdH;HJ
Exact match found for rate rule [Od_CO-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])C(=C)C[O](9813)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C(CO)C(=C)[O](28348)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]C(=CO)C(=C)O(28349)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C([O])C(=C)CO(28350)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1=C([O])CC1O(28351)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.6836e+08,'s^-1'), n=0.948854, Ea=(125.32,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1OCC=1[CH]O(28352)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C(O)C(=C)C=O(14411)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]C(=CO)C([CH2])[O](23305)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C([CH]O)=C(C)[O](23303)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]C([O])=C(C)[CH]O(28353)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(O)C(=C)[O](14497)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C1OC(O)C1=C(28293)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=[C]C(=C)[O](4767)', '[CH]O(5471)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C1([CH]O)CC1=O(28256)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(245.058,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(CO)_D;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(=[C]O)C(C)=O(28354)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['[CH2]C(=C[O])C(C)=O(14420)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(484628,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;XH_out] for rate rule [R5H_SSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(O)=C([CH2])[C]O(28355)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C(O)=C([CH2])[CH][O](28356)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=CO)C[C]=O(14236)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['O=C1CCC1=CO(28245)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['CH2(T)(28)', '[CH2]C([C]=O)=CO(28357)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]C(=CO)C(C)=O(28358)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C(C=O)C(C)=O(14422)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]C([CH]O)=C([CH2])O(23433)'],
    products = ['[CH2]C(=CO)C(=C)[O](14417)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C(=CO)C(=C)[O](14417)'],
    products = ['C=C1C(=O)CC1O(28255)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '5108',
    isomers = [
        '[CH2]C(=CO)C(=C)[O](14417)',
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
    label = '5108',
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


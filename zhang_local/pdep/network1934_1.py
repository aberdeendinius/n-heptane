species(
    label = '[CH2]C(C[O])CO[O](6836)',
    structure = SMILES('[CH2]C(C[O])CO[O]'),
    E0 = (181.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,1467.44,1468.01],'cm^-1')),
        HinderedRotor(inertia=(0.24356,'amu*angstrom^2'), symmetry=1, barrier=(5.59992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243682,'amu*angstrom^2'), symmetry=1, barrier=(5.60274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243564,'amu*angstrom^2'), symmetry=1, barrier=(5.60001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243846,'amu*angstrom^2'), symmetry=1, barrier=(5.6065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938385,0.0765282,-0.000117376,1.0704e-07,-3.79761e-11,21906.4,29.8831], Tmin=(100,'K'), Tmax=(867.715,'K')), NASAPolynomial(coeffs=[4.5136,0.0395006,-1.7849e-05,3.28422e-09,-2.20539e-13,22059.4,17.5992], Tmin=(867.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = 'C=O(355)',
    structure = SMILES('C=O'),
    E0 = (-118.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3229,-0.00506327,2.15156e-05,-1.76521e-08,4.31815e-12,-14279,2.39243], Tmin=(100,'K'), Tmax=(1402.28,'K')), NASAPolynomial(coeffs=[3.17995,0.00955599,-6.27302e-06,1.33555e-09,-9.68411e-14,-15075.2,4.31078], Tmin=(1402.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH)"""),
)

species(
    label = 'C=CCO[O](6082)',
    structure = SMILES('C=CCO[O]'),
    E0 = (76.4976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.316477,'amu*angstrom^2'), symmetry=1, barrier=(7.27644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315948,'amu*angstrom^2'), symmetry=1, barrier=(7.26427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3494.07,'J/mol'), sigma=(5.81539,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.77 K, Pc=40.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56676,0.0335761,-2.36362e-05,9.70939e-09,-1.82265e-12,9250.36,17.5875], Tmin=(100,'K'), Tmax=(1155.37,'K')), NASAPolynomial(coeffs=[5.58064,0.0231419,-1.00897e-05,1.8929e-09,-1.31328e-13,8553.92,2.61201], Tmin=(1155.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.4976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = 'C=C(C[O])CO[O](7331)',
    structure = SMILES('C=C(C[O])CO[O]'),
    E0 = (110.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,2479.76],'cm^-1')),
        HinderedRotor(inertia=(0.171294,'amu*angstrom^2'), symmetry=1, barrier=(3.93839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170163,'amu*angstrom^2'), symmetry=1, barrier=(3.91238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170817,'amu*angstrom^2'), symmetry=1, barrier=(3.92741,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35676,0.067573,-0.000105349,1.01585e-07,-3.82055e-11,13317.4,27.299], Tmin=(100,'K'), Tmax=(840.721,'K')), NASAPolynomial(coeffs=[2.65144,0.0404004,-1.9377e-05,3.68295e-09,-2.52769e-13,13842.3,25.694], Tmin=(840.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(C=O)CO[O](6958)',
    structure = SMILES('[CH2]C(C=O)CO[O]'),
    E0 = (41.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.29837,'amu*angstrom^2'), symmetry=1, barrier=(6.86011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298192,'amu*angstrom^2'), symmetry=1, barrier=(6.85601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298383,'amu*angstrom^2'), symmetry=1, barrier=(6.8604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298201,'amu*angstrom^2'), symmetry=1, barrier=(6.85623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751342,0.0797467,-0.000126333,1.12364e-07,-3.90682e-11,5079.74,27.7333], Tmin=(100,'K'), Tmax=(849.406,'K')), NASAPolynomial(coeffs=[7.09403,0.0335267,-1.58359e-05,2.9765e-09,-2.02518e-13,4592.1,1.64089], Tmin=(849.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(ROOJ)"""),
)

species(
    label = 'C=CC[O](423)',
    structure = SMILES('C=CC[O]'),
    E0 = (78.6931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,2140.08],'cm^-1')),
        HinderedRotor(inertia=(0.149731,'amu*angstrom^2'), symmetry=1, barrier=(3.44262,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31103,0.017642,3.66009e-07,-4.8614e-09,1.27726e-12,9487.15,14.2692], Tmin=(100,'K'), Tmax=(1683.2,'K')), NASAPolynomial(coeffs=[6.05036,0.0193257,-8.4361e-06,1.51683e-09,-9.96093e-14,7404.3,-3.82071], Tmin=(1683.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.6931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]O[O](61)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (200.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1055.29],'cm^-1')),
        HinderedRotor(inertia=(0.00752578,'amu*angstrom^2'), symmetry=1, barrier=(5.89392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3200.22,'J/mol'), sigma=(5.39124,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=499.87 K, Pc=46.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38991,0.0147959,-1.79906e-05,1.54065e-08,-5.43781e-12,24103.1,9.19559], Tmin=(100,'K'), Tmax=(865.677,'K')), NASAPolynomial(coeffs=[3.82258,0.00980114,-4.14553e-06,7.46963e-10,-4.99013e-14,24140.4,7.8189], Tmin=(865.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = '[CH2][O](357)',
    structure = SMILES('[CH2][O]'),
    E0 = (185.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81286,-0.00139668,2.72219e-05,-3.80869e-08,1.59436e-11,22322.5,7.76013], Tmin=(100,'K'), Tmax=(884.279,'K')), NASAPolynomial(coeffs=[6.98151,-0.0011443,2.05219e-06,-4.58191e-10,3.18315e-14,21191.8,-10.3616], Tmin=(884.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(CsJOH) + radical(H3COJ)"""),
)

species(
    label = '[CH2][CH]CO[O](6081)',
    structure = SMILES('[CH2][CH]CO[O]'),
    E0 = (348.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.00295793,'amu*angstrom^2'), symmetry=1, barrier=(5.15327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223054,'amu*angstrom^2'), symmetry=1, barrier=(5.12845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222928,'amu*angstrom^2'), symmetry=1, barrier=(5.12556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03778,0.0476625,-6.74371e-05,5.88574e-08,-2.07062e-11,42019.6,22.7181], Tmin=(100,'K'), Tmax=(834.77,'K')), NASAPolynomial(coeffs=[5.09101,0.024729,-1.13075e-05,2.11543e-09,-1.4439e-13,41799.2,10.2722], Tmin=(834.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'C[C](C[O])CO[O](7332)',
    structure = SMILES('C[C](C[O])CO[O]'),
    E0 = (167.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55232,0.0610603,-8.07941e-05,7.63519e-08,-2.97237e-11,20181.6,30.7406], Tmin=(100,'K'), Tmax=(807.833,'K')), NASAPolynomial(coeffs=[2.09679,0.0435128,-2.06347e-05,3.94746e-09,-2.73793e-13,20578.2,31.2294], Tmin=(807.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CJCOOH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([CH]O)CO[O](7333)',
    structure = SMILES('[CH2]C([CH]O)CO[O]'),
    E0 = (135.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330454,0.0877273,-0.000135217,1.12699e-07,-3.62995e-11,16469.2,30.7541], Tmin=(100,'K'), Tmax=(899.312,'K')), NASAPolynomial(coeffs=[9.74699,0.0300172,-1.25611e-05,2.20441e-09,-1.43008e-13,15415.5,-10.1176], Tmin=(899.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([CH]OO)C[O](7334)',
    structure = SMILES('[CH2]C([CH]OO)C[O]'),
    E0 = (217.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667436,0.0829894,-0.000129996,1.18465e-07,-4.21207e-11,26314.9,30.2459], Tmin=(100,'K'), Tmax=(854.774,'K')), NASAPolynomial(coeffs=[5.37683,0.0403046,-1.88588e-05,3.52681e-09,-2.39238e-13,26264.1,12.6768], Tmin=(854.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(CCsJOOH)"""),
)

species(
    label = 'CC([CH]O[O])C[O](7335)',
    structure = SMILES('CC([CH]O[O])C[O]'),
    E0 = (164.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.941751,0.0774103,-0.000121357,1.14622e-07,-4.22523e-11,19920.8,28.237], Tmin=(100,'K'), Tmax=(844.221,'K')), NASAPolynomial(coeffs=[3.60108,0.0426967,-2.03869e-05,3.86079e-09,-2.64223e-13,20259.8,20.5248], Tmin=(844.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'CC([CH][O])CO[O](7336)',
    structure = SMILES('CC([CH][O])CO[O]'),
    E0 = (156.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70555,0.0829638,-0.000134084,1.24619e-07,-4.48465e-11,18932.6,28.4853], Tmin=(100,'K'), Tmax=(853.623,'K')), NASAPolynomial(coeffs=[4.9269,0.0403862,-1.92068e-05,3.61664e-09,-2.46103e-13,19042.5,13.6527], Tmin=(853.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C](CO)CO[O](7337)',
    structure = SMILES('[CH2][C](CO)CO[O]'),
    E0 = (146.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18891,0.0656464,-8.10677e-05,6.28405e-08,-2.02083e-11,17717.7,32.9697], Tmin=(100,'K'), Tmax=(840.901,'K')), NASAPolynomial(coeffs=[6.99937,0.0330014,-1.3906e-05,2.51544e-09,-1.69046e-13,16917.5,6.99678], Tmin=(840.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(C2CJCOOH)"""),
)

species(
    label = '[CH2]C([CH]O[O])CO(7338)',
    structure = SMILES('[CH2]C([CH]O[O])CO'),
    E0 = (144.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56999,0.0821238,-0.000122251,1.02269e-07,-3.34475e-11,17457.3,30.4945], Tmin=(100,'K'), Tmax=(886.443,'K')), NASAPolynomial(coeffs=[8.44464,0.0322876,-1.37181e-05,2.44306e-09,-1.6067e-13,16623.1,-3.37734], Tmin=(886.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C[O])COO(7339)',
    structure = SMILES('[CH2][C](C[O])COO'),
    E0 = (220.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,292.948,978.026,2343.59],'cm^-1')),
        HinderedRotor(inertia=(0.0676406,'amu*angstrom^2'), symmetry=1, barrier=(3.33661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0676406,'amu*angstrom^2'), symmetry=1, barrier=(3.33661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0676406,'amu*angstrom^2'), symmetry=1, barrier=(3.33661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0676406,'amu*angstrom^2'), symmetry=1, barrier=(3.33661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0676406,'amu*angstrom^2'), symmetry=1, barrier=(3.33661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27947,0.0666139,-8.92931e-05,7.99051e-08,-2.94004e-11,26575.7,32.7446], Tmin=(100,'K'), Tmax=(817.769,'K')), NASAPolynomial(coeffs=[3.89228,0.0410861,-1.90863e-05,3.60861e-09,-2.48399e-13,26574.6,23.2711], Tmin=(817.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CJCOOH) + radical(CCOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH][O])COO(7340)',
    structure = SMILES('[CH2]C([CH][O])COO'),
    E0 = (209.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430544,0.0885538,-0.000142778,1.28567e-07,-4.47804e-11,25326.7,30.4964], Tmin=(100,'K'), Tmax=(864.71,'K')), NASAPolynomial(coeffs=[6.69643,0.0380048,-1.76849e-05,3.28415e-09,-2.21242e-13,25049.3,5.83962], Tmin=(864.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = '[CH2]C([CH2])C[O](454)',
    structure = SMILES('[CH2]C([CH2])C[O]'),
    E0 = (322.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1468.57,1474.61],'cm^-1')),
        HinderedRotor(inertia=(0.0980627,'amu*angstrom^2'), symmetry=1, barrier=(2.25466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0998327,'amu*angstrom^2'), symmetry=1, barrier=(2.29535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103728,'amu*angstrom^2'), symmetry=1, barrier=(2.38492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3512.21,'J/mol'), sigma=(6.17291,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.60 K, Pc=33.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07828,0.045281,-4.89163e-05,3.82982e-08,-1.26804e-11,38877.8,22.6574], Tmin=(100,'K'), Tmax=(892.843,'K')), NASAPolynomial(coeffs=[3.94866,0.0300676,-1.18763e-05,2.06845e-09,-1.35541e-13,38816.2,15.3713], Tmin=(892.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH2][CH]C[O](429)',
    structure = SMILES('[CH2][CH]C[O]'),
    E0 = (350.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2183.56,2183.57],'cm^-1')),
        HinderedRotor(inertia=(0.0017033,'amu*angstrom^2'), symmetry=1, barrier=(5.76299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301146,'amu*angstrom^2'), symmetry=1, barrier=(38.1405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29263,0.0206459,-6.09305e-06,-5.06452e-10,3.17593e-13,42176.5,17.7445], Tmin=(100,'K'), Tmax=(2029.72,'K')), NASAPolynomial(coeffs=[9.75428,0.0143081,-6.1363e-06,1.06034e-09,-6.66188e-14,38235.9,-21.2491], Tmin=(2029.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2][C](C[O])CO[O](7341)',
    structure = SMILES('[CH2][C](C[O])CO[O]'),
    E0 = (372.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,1775.75,1776.6],'cm^-1')),
        HinderedRotor(inertia=(0.141872,'amu*angstrom^2'), symmetry=1, barrier=(3.26191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141715,'amu*angstrom^2'), symmetry=1, barrier=(3.2583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145486,'amu*angstrom^2'), symmetry=1, barrier=(3.25869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141691,'amu*angstrom^2'), symmetry=1, barrier=(3.25776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51661,0.0633391,-9.50647e-05,9.04085e-08,-3.35132e-11,44847.2,32.6519], Tmin=(100,'K'), Tmax=(858.249,'K')), NASAPolynomial(coeffs=[2.37257,0.0393042,-1.80231e-05,3.35051e-09,-2.26809e-13,45438.6,32.9543], Tmin=(858.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(CCOJ) + radical(C2CJCOOH)"""),
)

species(
    label = '[CH2]C([CH]O[O])C[O](7342)',
    structure = SMILES('[CH2]C([CH]O[O])C[O]'),
    E0 = (369.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1264.86,1265.76],'cm^-1')),
        HinderedRotor(inertia=(0.178033,'amu*angstrom^2'), symmetry=1, barrier=(4.09332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17819,'amu*angstrom^2'), symmetry=1, barrier=(4.09693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178104,'amu*angstrom^2'), symmetry=1, barrier=(4.09496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178242,'amu*angstrom^2'), symmetry=1, barrier=(4.09813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921242,0.079515,-0.000135074,1.28129e-07,-4.59322e-11,44585.8,30.0937], Tmin=(100,'K'), Tmax=(879.901,'K')), NASAPolynomial(coeffs=[3.76086,0.038692,-1.78958e-05,3.2928e-09,-2.19672e-13,45166.7,22.898], Tmin=(879.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(Isobutyl) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH][O])CO[O](7343)',
    structure = SMILES('[CH2]C([CH][O])CO[O]'),
    E0 = (361.595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,180,1347.94,1349.19],'cm^-1')),
        HinderedRotor(inertia=(0.251838,'amu*angstrom^2'), symmetry=1, barrier=(5.79025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251858,'amu*angstrom^2'), symmetry=1, barrier=(5.79071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251687,'amu*angstrom^2'), symmetry=1, barrier=(5.78677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252044,'amu*angstrom^2'), symmetry=1, barrier=(5.79498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686803,0.0850524,-0.000147778,1.38166e-07,-4.85891e-11,43597.5,30.3354], Tmin=(100,'K'), Tmax=(887.727,'K')), NASAPolynomial(coeffs=[5.05704,0.0364329,-1.67457e-05,3.05583e-09,-2.02153e-13,43961.4,16.192], Tmin=(887.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=C(CO)CO[O](7344)',
    structure = SMILES('C=C(CO)CO[O]'),
    E0 = (-115.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0526,0.0695793,-9.01881e-05,7.23496e-08,-2.41212e-11,-13813.1,27.5338], Tmin=(100,'K'), Tmax=(801.797,'K')), NASAPolynomial(coeffs=[7.20436,0.0342298,-1.5339e-05,2.86708e-09,-1.96631e-13,-14649.8,0.148118], Tmin=(801.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CC(C=O)CO[O](7345)',
    structure = SMILES('CC(C=O)CO[O]'),
    E0 = (-169.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08821,0.0691008,-8.82889e-05,7.16583e-08,-2.45031e-11,-20247.8,26.1759], Tmin=(100,'K'), Tmax=(781.259,'K')), NASAPolynomial(coeffs=[6.61138,0.0363119,-1.66746e-05,3.1581e-09,-2.18501e-13,-20973.2,1.77415], Tmin=(781.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(C[O])COO(7346)',
    structure = SMILES('C=C(C[O])COO'),
    E0 = (-41.9928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1244,0.0707968,-9.94441e-05,9.10264e-08,-3.41538e-11,-4954.34,27.3743], Tmin=(100,'K'), Tmax=(798.917,'K')), NASAPolynomial(coeffs=[4.11999,0.0422728,-2.04938e-05,3.95399e-09,-2.7545e-13,-5001.34,16.2963], Tmin=(798.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.9928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(C=O)COO(7347)',
    structure = SMILES('[CH2]C(C=O)COO'),
    E0 = (-110.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533861,0.0827718,-0.000119608,1.00521e-07,-3.43474e-11,-13192.7,27.7566], Tmin=(100,'K'), Tmax=(793.707,'K')), NASAPolynomial(coeffs=[8.54442,0.0354327,-1.69733e-05,3.25262e-09,-2.25632e-13,-14244.8,-7.65615], Tmin=(793.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
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
    label = '[O]CC1COC1(7348)',
    structure = SMILES('[O]CC1COC1'),
    E0 = (-63.2635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30536,0.0324561,9.20836e-06,-3.20403e-08,1.52703e-11,-7543.21,18.3951], Tmin=(100,'K'), Tmax=(860.778,'K')), NASAPolynomial(coeffs=[6.28694,0.028374,-8.80675e-06,1.37421e-09,-8.68792e-14,-8762.88,-3.31989], Tmin=(860.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.2635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Oxetane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1COOC1(7349)',
    structure = SMILES('[CH2]C1COOC1'),
    E0 = (47.8401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11197,0.0266905,5.03128e-05,-9.45657e-08,4.32355e-11,5836.05,18.7892], Tmin=(100,'K'), Tmax=(865.857,'K')), NASAPolynomial(coeffs=[13.8373,0.0141064,7.46276e-08,-4.18584e-10,3.74686e-14,2246.78,-45.0914], Tmin=(865.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.8401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(Isobutyl)"""),
)

species(
    label = '[O]C[CH]CCO[O](7350)',
    structure = SMILES('[O]C[CH]CCO[O]'),
    E0 = (185.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,1817.92,1821.66],'cm^-1')),
        HinderedRotor(inertia=(0.189554,'amu*angstrom^2'), symmetry=1, barrier=(4.35821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18949,'amu*angstrom^2'), symmetry=1, barrier=(4.35676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187729,'amu*angstrom^2'), symmetry=1, barrier=(4.31625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18946,'amu*angstrom^2'), symmetry=1, barrier=(4.35605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09998,0.0722243,-0.000106908,9.88871e-08,-3.64432e-11,22365.6,30.2962], Tmin=(100,'K'), Tmax=(830.131,'K')), NASAPolynomial(coeffs=[3.99519,0.0410594,-1.94892e-05,3.70257e-09,-2.54699e-13,22478.1,20.44], Tmin=(830.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCJCO) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]CC[CH]CO[O](6835)',
    structure = SMILES('[O]CC[CH]CO[O]'),
    E0 = (185.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,1796.6,1796.82],'cm^-1')),
        HinderedRotor(inertia=(0.236016,'amu*angstrom^2'), symmetry=1, barrier=(5.42647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235468,'amu*angstrom^2'), symmetry=1, barrier=(5.41387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23562,'amu*angstrom^2'), symmetry=1, barrier=(5.41736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23612,'amu*angstrom^2'), symmetry=1, barrier=(5.42887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922927,0.0795637,-0.000131823,1.28072e-07,-4.77284e-11,22430.4,30.7395], Tmin=(100,'K'), Tmax=(850.4,'K')), NASAPolynomial(coeffs=[2.74718,0.0439458,-2.13073e-05,4.04714e-09,-2.76765e-13,23097.7,27.9822], Tmin=(850.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCJCOOH) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = '[CH2]C(=C)C[O](6285)',
    structure = SMILES('[CH2]C(=C)C[O]'),
    E0 = (191.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,647.829,649.857],'cm^-1')),
        HinderedRotor(inertia=(0.160033,'amu*angstrom^2'), symmetry=1, barrier=(3.67947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0123324,'amu*angstrom^2'), symmetry=1, barrier=(3.69332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.65316,0.0333255,-1.65856e-05,3.65186e-09,-3.087e-13,23032.8,16.9441], Tmin=(100,'K'), Tmax=(2703.22,'K')), NASAPolynomial(coeffs=[16.127,0.0133877,-5.52202e-06,9.23318e-10,-5.63535e-14,15748.4,-61.4582], Tmin=(2703.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = '[O]OCC1COC1(6840)',
    structure = SMILES('[O]OCC1COC1'),
    E0 = (-65.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825903,0.0568855,-4.38266e-05,1.94992e-08,-3.45317e-12,-7747.56,24.3641], Tmin=(100,'K'), Tmax=(1570.89,'K')), NASAPolynomial(coeffs=[11.0824,0.0240685,-6.09209e-06,7.69715e-10,-4.0323e-14,-10143.1,-27.1181], Tmin=(1570.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Oxetane) + radical(ROOJ)"""),
)

species(
    label = '[O]CC1COOC1(7351)',
    structure = SMILES('[O]CC1COOC1'),
    E0 = (-91.3645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70264,0.0424211,3.16374e-06,-3.50708e-08,1.79438e-11,-10898.3,22.0357], Tmin=(100,'K'), Tmax=(892.765,'K')), NASAPolynomial(coeffs=[10.2384,0.0265604,-7.79565e-06,1.19673e-09,-7.63963e-14,-13314.4,-23.1722], Tmin=(892.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.3645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1COOOC1(7352)',
    structure = SMILES('[CH2]C1COOOC1'),
    E0 = (76.6331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28455,0.0468678,1.15339e-05,-6.19694e-08,3.3171e-11,9327.02,19.0285], Tmin=(100,'K'), Tmax=(859.286,'K')), NASAPolynomial(coeffs=[16.3585,0.0152674,-6.31234e-07,-2.95667e-10,3.03568e-14,5312.52,-59.6942], Tmin=(859.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(123trioxane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C[O])C[O](7353)',
    structure = SMILES('[CH2]C(C[O])C[O]'),
    E0 = (183.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1663.84,1663.84,1664.01],'cm^-1')),
        HinderedRotor(inertia=(0.166882,'amu*angstrom^2'), symmetry=1, barrier=(3.83694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166745,'amu*angstrom^2'), symmetry=1, barrier=(3.83378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166652,'amu*angstrom^2'), symmetry=1, barrier=(3.83166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71867,0.0603811,-9.36422e-05,9.43155e-08,-3.63279e-11,22141.3,25.7284], Tmin=(100,'K'), Tmax=(858.962,'K')), NASAPolynomial(coeffs=[0.22038,0.0427507,-1.98826e-05,3.71658e-09,-2.52171e-13,23306.5,38.0133], Tmin=(858.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]C[CH]CO[O](6847)',
    structure = SMILES('[O]C[CH]CO[O]'),
    E0 = (209.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,2503.52],'cm^-1')),
        HinderedRotor(inertia=(0.125425,'amu*angstrom^2'), symmetry=1, barrier=(2.88376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127123,'amu*angstrom^2'), symmetry=1, barrier=(2.92282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125858,'amu*angstrom^2'), symmetry=1, barrier=(2.89373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61109,0.0642626,-0.000115985,1.17039e-07,-4.40917e-11,25265.9,26.0252], Tmin=(100,'K'), Tmax=(862.72,'K')), NASAPolynomial(coeffs=[1.72127,0.0357136,-1.75976e-05,3.33966e-09,-2.27198e-13,26290.3,31.5572], Tmin=(862.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C([CH2])CO[O](6921)',
    structure = SMILES('[CH2]C([CH2])CO[O]'),
    E0 = (320.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.0985179,'amu*angstrom^2'), symmetry=1, barrier=(2.26512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985197,'amu*angstrom^2'), symmetry=1, barrier=(2.26516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0984856,'amu*angstrom^2'), symmetry=1, barrier=(2.26438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986061,'amu*angstrom^2'), symmetry=1, barrier=(2.26715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26784,0.061818,-7.41934e-05,5.33286e-08,-1.54728e-11,38644.1,26.2251], Tmin=(100,'K'), Tmax=(935.888,'K')), NASAPolynomial(coeffs=[8.27539,0.0267562,-9.80549e-06,1.627e-09,-1.03136e-13,37556.3,-5.92207], Tmin=(935.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(C[O])CO[O](7354)',
    structure = SMILES('[CH]C(C[O])CO[O]'),
    E0 = (424.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889882,0.077669,-0.000123788,1.13876e-07,-4.08683e-11,51150.1,29.2784], Tmin=(100,'K'), Tmax=(845.992,'K')), NASAPolynomial(coeffs=[5.34658,0.0372598,-1.78534e-05,3.37839e-09,-2.30848e-13,51088.1,12.6125], Tmin=(845.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(ROOJ) + radical(CCJ2_triplet)"""),
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
    E0 = (181.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (321.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (270.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (315.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (275.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (289.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (283.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (270.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (347.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (298.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (299.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (256.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (290.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (304.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (251.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (325.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (550.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (534.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (584.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (581.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (573.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (244.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (244.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (206.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (206.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (264.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (290.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (430.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (431.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (320.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (189.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (188.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (190.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (431.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (625.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (727.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (636.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['C=O(355)', 'C=CCO[O](6082)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=C(C[O])CO[O](7331)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.00507518,'m^3/(mol*s)'), n=2.82235, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C(C=O)CO[O](6958)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC[O](423)', '[CH2]O[O](61)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.00140605,'m^3/(mol*s)'), n=2.41796, Ea=(36.1629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ] for rate rule [Cds-CsH_Cds-HH;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][O](357)', 'C=CCO[O](6082)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.540096,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=O(355)', '[CH2][CH]CO[O](6081)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(225.36,'m^3/(mol*s)'), n=0.996465, Ea=(58.8821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['C[C](C[O])CO[O](7332)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2]C([CH]O)CO[O](7333)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2]C([CH]OO)C[O](7334)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['CC([CH]O[O])C[O](7335)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['CC([CH][O])CO[O](7336)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2][C](CO)CO[O](7337)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2]C([CH]O[O])CO(7338)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.22e+08,'s^-1'), n=1.09, Ea=(109.37,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](C[O])COO(7339)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH][O])COO(7340)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O2(2)', '[CH2]C([CH2])C[O](454)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.66237e+08,'m^3/(mol*s)'), n=-0.783071, Ea=(11.596,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for O2_birad;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C[O](429)', '[CH2]O[O](61)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.88087e+06,'m^3/(mol*s)'), n=0.114385, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_pri_rad] for rate rule [Y_rad;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][O](357)', '[CH2][CH]CO[O](6081)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2][C](C[O])CO[O](7341)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2]C([CH]O[O])C[O](7342)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C([CH][O])CO[O](7343)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['C=C(CO)CO[O](7344)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['CC(C=O)CO[O](7345)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['C=C(C[O])COO(7346)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2]C(C=O)COO(7347)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['O(T)(63)', '[O]CC1COC1(7348)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['O(T)(63)', '[CH2]C1COOC1(7349)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(109.578,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;Y_rad_intra;OO] for rate rule [R4OO_SSS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation
Ea raised from 109.4 to 109.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C[CH]CCO[O](7350)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]CC[CH]CO[O](6835)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[O]O(16)', '[CH2]C(=C)C[O](6285)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[O]OCC1COC1(6840)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[O]CC1COOC1(7351)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C[O])CO[O](6836)'],
    products = ['[CH2]C1COOOC1(7352)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(T)(63)', '[CH2]C(C[O])C[O](7353)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(T)(28)', '[O]C[CH]CO[O](6847)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH2]C([CH2])CO[O](6921)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(187219,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', '[CH]C(C[O])CO[O](7354)'],
    products = ['[CH2]C(C[O])CO[O](6836)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '1934',
    isomers = [
        '[CH2]C(C[O])CO[O](6836)',
    ],
    reactants = [
        ('C=O(355)', 'C=CCO[O](6082)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '1934',
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


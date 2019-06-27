species(
    label = '[CH]=C=COC([O])O[O](22517)',
    structure = SMILES('[CH]=C=COC([O])O[O]'),
    E0 = (216.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,280.085],'cm^-1')),
        HinderedRotor(inertia=(0.105986,'amu*angstrom^2'), symmetry=1, barrier=(2.43684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471398,'amu*angstrom^2'), symmetry=1, barrier=(10.8384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.34189,'amu*angstrom^2'), symmetry=1, barrier=(53.8447,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725813,0.0823308,-0.000143902,1.33856e-07,-4.7845e-11,26125.6,27.791], Tmin=(100,'K'), Tmax=(842.231,'K')), NASAPolynomial(coeffs=[6.90433,0.0318694,-1.64209e-05,3.17807e-09,-2.19033e-13,25833.9,3.49054], Tmin=(842.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(OCOJ) + radical(ROOJ) + radical(C=C=CJ)"""),
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
    label = '[CH]=[C]C1OC(O[O])O1(24712)',
    structure = SMILES('[CH]=[C]C1OC(O[O])O1'),
    E0 = (342.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43838,0.0633388,-9.56008e-05,9.10231e-08,-3.60175e-11,41259.6,26.5067], Tmin=(100,'K'), Tmax=(734.514,'K')), NASAPolynomial(coeffs=[4.74635,0.0356669,-1.93683e-05,3.93194e-09,-2.82455e-13,41034.2,13.3417], Tmin=(734.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = '[CH]=[C]C1OOC([O])O1(24713)',
    structure = SMILES('[CH]=[C]C1OOC([O])O1'),
    E0 = (328.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97031,0.0467971,-4.05855e-05,2.23729e-08,-5.53944e-12,39584,24.8556], Tmin=(100,'K'), Tmax=(931.801,'K')), NASAPolynomial(coeffs=[6.02627,0.0293861,-1.2558e-05,2.32066e-09,-1.59565e-13,38828.1,5.57441], Tmin=(931.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(124trioxolane) + radical(Cds_S) + radical(OCOJ) + radical(Cds_P)"""),
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
    label = '[CH]=C=COC(=O)O[O](24714)',
    structure = SMILES('[CH]=C=COC(=O)O[O]'),
    E0 = (-44.2995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,3010,987.5,1337.5,450,1655,540,610,2055,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5708,'amu*angstrom^2'), symmetry=1, barrier=(36.1158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57421,'amu*angstrom^2'), symmetry=1, barrier=(36.1941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57199,'amu*angstrom^2'), symmetry=1, barrier=(36.1432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.943503,0.0783049,-8.77e-05,4.32936e-08,-7.68056e-12,-5123.71,32.297], Tmin=(100,'K'), Tmax=(1668.37,'K')), NASAPolynomial(coeffs=[24.4928,-0.00340582,4.39885e-06,-9.54376e-10,6.56426e-14,-10726.6,-94.7932], Tmin=(1668.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.2995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C(=O)OOJ) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=COC=O(18779)',
    structure = SMILES('[CH]=C=COC=O'),
    E0 = (0.856013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.53408,'amu*angstrom^2'), symmetry=1, barrier=(35.2716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53438,'amu*angstrom^2'), symmetry=1, barrier=(35.2784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65569,0.0438498,-2.93891e-05,-1.14417e-09,5.77321e-12,194.379,19.0267], Tmin=(100,'K'), Tmax=(950.188,'K')), NASAPolynomial(coeffs=[14.5183,0.0082171,-2.36647e-06,4.03443e-10,-2.95137e-14,-3085.82,-46.7689], Tmin=(950.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.856013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-OdOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CO[C](O)O[O](24715)',
    structure = SMILES('[CH]=C=CO[C](O)O[O]'),
    E0 = (194.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220217,0.0857196,-0.000128955,9.50384e-08,-2.70203e-11,23487.2,29.8233], Tmin=(100,'K'), Tmax=(871.385,'K')), NASAPolynomial(coeffs=[15.7617,0.0143839,-6.16868e-06,1.1062e-09,-7.33898e-14,20778.4,-43.0171], Tmin=(871.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH][C]=COC(=O)OO(24716)',
    structure = SMILES('[CH][C]=COC(=O)OO'),
    E0 = (38.7571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0482216,0.070363,-3.81331e-05,-1.7441e-08,1.57742e-11,4818.7,29.8537], Tmin=(100,'K'), Tmax=(959.57,'K')), NASAPolynomial(coeffs=[24.4524,0.00803805,-2.30422e-06,4.61869e-10,-3.97005e-14,-1678.95,-96.3282], Tmin=(959.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.7571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C#COC([O])O[O](24717)',
    structure = SMILES('[CH2]C#COC([O])O[O]'),
    E0 = (244.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2100,2250,500,550,3000,3100,440,815,1455,1000,180,180,180,1349.83],'cm^-1')),
        HinderedRotor(inertia=(0.349688,'amu*angstrom^2'), symmetry=1, barrier=(8.04003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312336,'amu*angstrom^2'), symmetry=1, barrier=(40.4864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3502,'amu*angstrom^2'), symmetry=1, barrier=(8.0518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312779,'amu*angstrom^2'), symmetry=1, barrier=(40.506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996911,0.0790306,-0.000146336,1.48069e-07,-5.70405e-11,29540.6,27.2609], Tmin=(100,'K'), Tmax=(821.036,'K')), NASAPolynomial(coeffs=[3.19401,0.0396024,-2.18241e-05,4.3563e-09,-3.06194e-13,30148,22.9903], Tmin=(821.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Propargyl) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C=[C]OC(O)O[O](24718)',
    structure = SMILES('[CH]=C=[C]OC(O)O[O]'),
    E0 = (228.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,1685,370,184.694,189.522],'cm^-1')),
        HinderedRotor(inertia=(0.319732,'amu*angstrom^2'), symmetry=1, barrier=(8.62801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338161,'amu*angstrom^2'), symmetry=1, barrier=(8.62786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32735,'amu*angstrom^2'), symmetry=1, barrier=(8.67265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878806,'amu*angstrom^2'), symmetry=1, barrier=(22.8989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360064,0.0871742,-0.000147121,1.24205e-07,-4.01286e-11,27626.9,31.7592], Tmin=(100,'K'), Tmax=(877.8,'K')), NASAPolynomial(coeffs=[11.8989,0.0206766,-9.70782e-06,1.78207e-09,-1.18163e-13,26137.3,-19.3511], Tmin=(877.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=[C]OC([O])OO(24719)',
    structure = SMILES('[CH]=C=[C]OC([O])OO'),
    E0 = (304.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1310,387.5,850,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,180,180,2756.94],'cm^-1')),
        HinderedRotor(inertia=(2.04303,'amu*angstrom^2'), symmetry=1, barrier=(46.9732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04266,'amu*angstrom^2'), symmetry=1, barrier=(46.9648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04259,'amu*angstrom^2'), symmetry=1, barrier=(46.9631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132848,'amu*angstrom^2'), symmetry=1, barrier=(3.05443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660189,0.0878426,-0.000166804,1.64587e-07,-6.10304e-11,36676.7,30.3947], Tmin=(100,'K'), Tmax=(847.081,'K')), NASAPolynomial(coeffs=[4.45417,0.0374763,-2.01526e-05,3.94494e-09,-2.72732e-13,37198.2,19.5926], Tmin=(847.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(OCOJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C][CH]OC(=O)O[O](24720)',
    structure = SMILES('C=[C][CH]OC(=O)O[O]'),
    E0 = (5.56645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09095,0.0541715,-3.54873e-05,2.94452e-09,3.45495e-12,782.696,30.3953], Tmin=(100,'K'), Tmax=(1052.43,'K')), NASAPolynomial(coeffs=[16.0448,0.0156855,-6.78753e-06,1.33148e-09,-9.72987e-14,-3381.08,-47.3404], Tmin=(1052.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.56645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(C(=O)OOJ) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
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
    label = '[O]OC([O])[O](4172)',
    structure = SMILES('[O]OC([O])[O]'),
    E0 = (77.6018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2220.26,2220.69],'cm^-1')),
        HinderedRotor(inertia=(3.41843e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.0162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3674,0.0423788,-0.000145172,1.9796e-07,-8.61253e-11,9328.67,15.7326], Tmin=(100,'K'), Tmax=(861.03,'K')), NASAPolynomial(coeffs=[-15.5083,0.0534406,-3.09507e-05,6.16376e-09,-4.27232e-13,15419.6,120.467], Tmin=(861.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.6018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(ROOJ) + radical(OCOJ)"""),
)

species(
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH][C]=COC=O(19586)',
    structure = SMILES('[CH][C]=COC=O'),
    E0 = (278.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1685,370,200.4,200.401,200.403,200.411,200.412],'cm^-1')),
        HinderedRotor(inertia=(1.72584,'amu*angstrom^2'), symmetry=1, barrier=(49.1886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72591,'amu*angstrom^2'), symmetry=1, barrier=(49.1887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72589,'amu*angstrom^2'), symmetry=1, barrier=(49.1886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45492,0.0489258,-4.01566e-05,1.63651e-08,-2.65594e-12,33569.4,21.6602], Tmin=(100,'K'), Tmax=(1471.71,'K')), NASAPolynomial(coeffs=[13.531,0.0161044,-6.7046e-06,1.21194e-09,-8.19001e-14,30014.9,-41.2664], Tmin=(1471.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=COC(=O)O[O](24721)',
    structure = SMILES('[CH][C]=COC(=O)O[O]'),
    E0 = (233.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336961,0.0664852,-4.18783e-05,-9.30577e-09,1.24799e-11,28185.7,29.576], Tmin=(100,'K'), Tmax=(949.266,'K')), NASAPolynomial(coeffs=[22.5945,0.00684038,-1.58198e-06,2.84988e-10,-2.48908e-14,22421.7,-84.7477], Tmin=(949.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(C(=O)OOJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]OC([O])O[O](24722)',
    structure = SMILES('[CH]=C=[C]OC([O])O[O]'),
    E0 = (456.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.15389,'amu*angstrom^2'), symmetry=1, barrier=(3.53824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153219,'amu*angstrom^2'), symmetry=1, barrier=(3.5228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152966,'amu*angstrom^2'), symmetry=1, barrier=(3.517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921786,0.0842617,-0.000171428,1.7351e-07,-6.44414e-11,54947.2,30.2155], Tmin=(100,'K'), Tmax=(866.432,'K')), NASAPolynomial(coeffs=[2.84451,0.035854,-1.91844e-05,3.70975e-09,-2.53073e-13,56097.9,29.7779], Tmin=(866.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CJO) + radical(C=C=CJ) + radical(OCOJ)"""),
)

species(
    label = '[O]OC([O])OC1[C]=C1(24723)',
    structure = SMILES('[O]OC([O])OC1[C]=C1'),
    E0 = (415.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.8699,0.0297431,-1.16879e-05,8.11554e-11,2.75974e-13,49773.1,11.5161], Tmin=(100,'K'), Tmax=(2775.01,'K')), NASAPolynomial(coeffs=[68.8864,-0.0342413,8.38469e-06,-1.25437e-09,8.21734e-14,4460.75,-375.447], Tmin=(2775.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-OsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(OCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]C1=COC(O[O])O1(24574)',
    structure = SMILES('[CH]C1=COC(O[O])O1'),
    E0 = (78.9056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1111,0.0405786,3.86326e-05,-9.38584e-08,4.27014e-11,9615.28,23.4546], Tmin=(100,'K'), Tmax=(937.408,'K')), NASAPolynomial(coeffs=[22.7231,0.00659716,-1.83322e-07,2.26764e-11,-1.13609e-14,3004.59,-93.0624], Tmin=(937.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.9056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
)

species(
    label = '[CH]C1=COC([O])OO1(24724)',
    structure = SMILES('[CH]C1=COC([O])OO1'),
    E0 = (190.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76882,0.047558,-2.80115e-05,7.83995e-09,-8.89417e-13,22985,22.2383], Tmin=(100,'K'), Tmax=(1934.5,'K')), NASAPolynomial(coeffs=[12.4224,0.0255294,-1.09306e-05,1.9535e-09,-1.28694e-13,18863.2,-36.1885], Tmin=(1934.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(124trioxene) + radical(OCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=COC(=O)OO(22526)',
    structure = SMILES('[CH]=C=COC(=O)OO'),
    E0 = (-238.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.128495,0.06653,-3.08168e-05,-3.18715e-08,2.35387e-11,-28550.9,27.6637], Tmin=(100,'K'), Tmax=(926.966,'K')), NASAPolynomial(coeffs=[27.3643,-0.00290061,3.70647e-06,-7.27151e-10,4.33614e-14,-35666.6,-112.813], Tmin=(926.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = 'C=C=COC(=O)O[O](24725)',
    structure = SMILES('C=C=COC(=O)O[O]'),
    E0 = (-198.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406843,0.0616317,-2.66058e-05,-2.94847e-08,2.08434e-11,-23761.9,26.7965], Tmin=(100,'K'), Tmax=(940.876,'K')), NASAPolynomial(coeffs=[25.0668,0.000353588,1.6412e-06,-2.92655e-10,1.20842e-14,-30330.4,-100.917], Tmin=(940.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C(=O)OOJ)"""),
)

species(
    label = '[CH]=[C]CO[C]([O])O[O](24726)',
    structure = SMILES('[CH]=[C]CO[C]([O])O[O]'),
    E0 = (633.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,360,370,350,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,180,1842.25,1844.74,1847.18],'cm^-1')),
        HinderedRotor(inertia=(0.217106,'amu*angstrom^2'), symmetry=1, barrier=(4.9917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217196,'amu*angstrom^2'), symmetry=1, barrier=(4.99377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218375,'amu*angstrom^2'), symmetry=1, barrier=(5.02086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21641,'amu*angstrom^2'), symmetry=1, barrier=(4.97569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06919,0.0920509,-0.000216218,2.42922e-07,-9.61477e-11,76228.8,32.0359], Tmin=(100,'K'), Tmax=(863.025,'K')), NASAPolynomial(coeffs=[-5.98494,0.0550364,-3.07241e-05,6.03848e-09,-4.15573e-13,80042.4,80.0691], Tmin=(863.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(ROOJ) + radical(Cds_S) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = '[CH][CH][CH]OC(=O)O[O](24727)',
    structure = SMILES('[CH][CH][CH]OC(=O)O[O]'),
    E0 = (365.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3050,390,425,1340,1360,335,370,180,1047.69,1047.69,1047.69,1047.69,1047.69,1047.69,1047.69,1047.69,1047.69,2259.49],'cm^-1')),
        HinderedRotor(inertia=(0.1086,'amu*angstrom^2'), symmetry=1, barrier=(2.49693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1086,'amu*angstrom^2'), symmetry=1, barrier=(2.49693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1086,'amu*angstrom^2'), symmetry=1, barrier=(2.49693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1086,'amu*angstrom^2'), symmetry=1, barrier=(2.49693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1086,'amu*angstrom^2'), symmetry=1, barrier=(2.49693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315562,0.0679246,-5.19549e-05,1.71922e-09,8.55437e-12,44103,33.8826], Tmin=(100,'K'), Tmax=(954.306,'K')), NASAPolynomial(coeffs=[23.0282,0.00397075,-5.44997e-07,1.15466e-10,-1.38425e-14,38345.2,-82.085], Tmin=(954.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsOs) + radical(CCsJOC(O)) + radical(CCJCO) + radical(C(=O)OOJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH]=C=COC1OO1(24728)',
    structure = SMILES('[CH]=C=COC1OO1'),
    E0 = (159.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2175,0.0467762,-8.09002e-06,-4.25791e-08,2.59723e-11,19312.5,21.584], Tmin=(100,'K'), Tmax=(891.481,'K')), NASAPolynomial(coeffs=[20.8213,-0.000840832,4.1487e-06,-9.68472e-10,6.77501e-14,14214.1,-79.7322], Tmin=(891.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(dioxirane) + radical(C=C=CJ)"""),
)

species(
    label = '[O]C1OC=C=CO1(24729)',
    structure = SMILES('[O]C1OC=C=CO1'),
    E0 = (-102.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.323257,0.101573,-0.000147445,1.06216e-07,-3.01489e-11,-12149,4.48398], Tmin=(100,'K'), Tmax=(864.037,'K')), NASAPolynomial(coeffs=[16.3735,0.0242776,-1.32578e-05,2.68246e-09,-1.92718e-13,-15034.3,-73.6281], Tmin=(864.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(OCOJ)"""),
)

species(
    label = '[CH]=C=COC1OOO1(24730)',
    structure = SMILES('[CH]=C=COC1OOO1'),
    E0 = (200.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1169,0.0429298,1.96579e-05,-7.0493e-08,3.34725e-11,24277.9,25.0376], Tmin=(100,'K'), Tmax=(949.832,'K')), NASAPolynomial(coeffs=[22.4033,0.00439646,-2.02919e-07,9.80928e-11,-1.82026e-14,17928.7,-88.6978], Tmin=(949.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclobutane) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C(C=O)C([O])O[O](22519)',
    structure = SMILES('[CH]=C(C=O)C([O])O[O]'),
    E0 = (211.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.265418,'amu*angstrom^2'), symmetry=1, barrier=(6.10249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266002,'amu*angstrom^2'), symmetry=1, barrier=(6.11592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264975,'amu*angstrom^2'), symmetry=1, barrier=(6.09229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4576.66,'J/mol'), sigma=(7.08346,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=714.86 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80174,0.0845562,-0.000162091,1.63317e-07,-6.20966e-11,25547.4,29.9109], Tmin=(100,'K'), Tmax=(829.38,'K')), NASAPolynomial(coeffs=[3.81262,0.03841,-2.14355e-05,4.28111e-09,-3.00191e-13,26135.7,22.5057], Tmin=(829.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=COC([O])[O](24731)',
    structure = SMILES('[CH]=C=COC([O])[O]'),
    E0 = (220.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,435.483,435.979,437.23,439.547],'cm^-1')),
        HinderedRotor(inertia=(0.000892829,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000895434,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.34608,0.0270737,-9.96393e-06,-4.64661e-11,2.55141e-13,26378.2,9.1375], Tmin=(100,'K'), Tmax=(2744.81,'K')), NASAPolynomial(coeffs=[55.0083,-0.0230408,5.25922e-06,-7.89544e-10,5.37348e-14,-9269.06,-295.872], Tmin=(2744.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[CH]=C=CO[CH]O[O](19776)',
    structure = SMILES('[CH]=C=CO[CH]O[O]'),
    E0 = (367.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18002,'amu*angstrom^2'), symmetry=1, barrier=(27.1309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1774,'amu*angstrom^2'), symmetry=1, barrier=(27.0707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17911,'amu*angstrom^2'), symmetry=1, barrier=(27.11,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.132685,0.0801725,-0.000106213,6.37141e-08,-1.40064e-11,44371.4,26.3834], Tmin=(100,'K'), Tmax=(1282.6,'K')), NASAPolynomial(coeffs=[21.8132,-0.000845586,3.2456e-06,-8.24662e-10,6.33117e-14,39776.2,-80.9225], Tmin=(1282.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[C]=C=COC([O])O[O](24732)',
    structure = SMILES('[C]=C=COC([O])O[O]'),
    E0 = (619.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.140706,'amu*angstrom^2'), symmetry=1, barrier=(3.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141835,'amu*angstrom^2'), symmetry=1, barrier=(3.26106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76406,'amu*angstrom^2'), symmetry=1, barrier=(40.5592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688592,0.0876153,-0.000173153,1.72077e-07,-6.39641e-11,74667.6,27.6733], Tmin=(100,'K'), Tmax=(843.272,'K')), NASAPolynomial(coeffs=[4.94592,0.034273,-1.93056e-05,3.83524e-09,-2.66899e-13,75128.1,14.848], Tmin=(843.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CdCdJ2_triplet) + radical(ROOJ) + radical(OCOJ)"""),
)

species(
    label = 'C#CCO[C]([O])O[O](24733)',
    structure = SMILES('C#CCO[C]([O])O[O]'),
    E0 = (314.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,750,770,3400,2100,2175,525,2750,2850,1437.5,1250,1305,750,350,180,1262.34,1262.41,1262.58],'cm^-1')),
        HinderedRotor(inertia=(0.208933,'amu*angstrom^2'), symmetry=1, barrier=(4.80378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209366,'amu*angstrom^2'), symmetry=1, barrier=(4.81374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20872,'amu*angstrom^2'), symmetry=1, barrier=(4.79888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209064,'amu*angstrom^2'), symmetry=1, barrier=(4.8068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13679,0.0869295,-0.000193891,2.12755e-07,-8.29151e-11,37864.5,28.5358], Tmin=(100,'K'), Tmax=(869.333,'K')), NASAPolynomial(coeffs=[-3.98004,0.0503473,-2.70253e-05,5.23136e-09,-3.56899e-13,41026.2,65.5723], Tmin=(869.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(OCOJ) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[C]#CCOC([O])O[O](24734)',
    structure = SMILES('[C]#CCOC([O])O[O]'),
    E0 = (446.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2175,525,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12464,0.0896913,-0.000205712,2.26723e-07,-8.78237e-11,53726.3,28.6348], Tmin=(100,'K'), Tmax=(879.873,'K')), NASAPolynomial(coeffs=[-5.27328,0.0519221,-2.73502e-05,5.22414e-09,-3.52385e-13,57440,73.3883], Tmin=(879.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(Acetyl) + radical(OCOJ)"""),
)

species(
    label = '[C]=C=COC(O)O[O](24735)',
    structure = SMILES('[C]=C=COC(O)O[O]'),
    E0 = (392.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.630295,'amu*angstrom^2'), symmetry=1, barrier=(14.4917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63042,'amu*angstrom^2'), symmetry=1, barrier=(14.4946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.630026,'amu*angstrom^2'), symmetry=1, barrier=(14.4855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.630055,'amu*angstrom^2'), symmetry=1, barrier=(14.4862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.181713,0.0898203,-0.000146068,1.18692e-07,-3.76795e-11,47344.9,29.0241], Tmin=(100,'K'), Tmax=(792.37,'K')), NASAPolynomial(coeffs=[13.8636,0.0193413,-9.97619e-06,1.94331e-09,-1.35016e-13,45221,-33.5193], Tmin=(792.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[C]=C=COC([O])OO(24736)',
    structure = SMILES('[C]=C=COC([O])OO'),
    E0 = (467.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,180,1638.18],'cm^-1')),
        HinderedRotor(inertia=(1.43385,'amu*angstrom^2'), symmetry=1, barrier=(32.9669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44151,'amu*angstrom^2'), symmetry=1, barrier=(33.1432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30436,'amu*angstrom^2'), symmetry=1, barrier=(6.99784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43939,'amu*angstrom^2'), symmetry=1, barrier=(33.0943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435673,0.0911023,-0.000168268,1.62982e-07,-6.05966e-11,56396.7,27.821], Tmin=(100,'K'), Tmax=(817.001,'K')), NASAPolynomial(coeffs=[6.47539,0.0360361,-2.03568e-05,4.09035e-09,-2.88231e-13,56260.7,5.11119], Tmin=(817.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CdCdJ2_triplet) + radical(OCOJ)"""),
)

species(
    label = '[O]OC1O[CH][C]=CO1(24608)',
    structure = SMILES('[O]OC1O[CH][C]=CO1'),
    E0 = (106.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53287,0.0454135,-1.83722e-05,-1.35357e-08,1.00045e-11,12911.9,20.8924], Tmin=(100,'K'), Tmax=(948.997,'K')), NASAPolynomial(coeffs=[13.733,0.0154466,-4.92017e-06,8.38873e-10,-5.85021e-14,9630.2,-42.4179], Tmin=(948.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(24dihydro13dioxin) + radical(Cds_S) + radical(C=CCJ(O)C) + radical(ROOJ)"""),
)

species(
    label = '[O]C1OC=[C][CH]OO1(24737)',
    structure = SMILES('[O]C1OC=[C][CH]OO1'),
    E0 = (185.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7848,0.000114392,0.000184831,-2.74811e-07,1.14868e-10,22395.7,23.8369], Tmin=(100,'K'), Tmax=(918.099,'K')), NASAPolynomial(coeffs=[35.9486,-0.0227764,1.6443e-05,-3.10871e-09,1.92961e-13,10814.2,-166.974], Tmin=(918.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(OCOJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = 'C#CCOC(=O)O[O](24738)',
    structure = SMILES('C#CCOC(=O)O[O]'),
    E0 = (-177.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943312,0.0542095,-2.29007e-05,-2.19859e-08,1.59211e-11,-21198,27.13], Tmin=(100,'K'), Tmax=(935.845,'K')), NASAPolynomial(coeffs=[19.423,0.0080343,-1.47986e-06,2.1817e-10,-1.83936e-14,-26093.6,-68.4749], Tmin=(935.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Cds-OdOsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C(=O)OOJ)"""),
)

species(
    label = 'C#CC1OC([O])O1(24739)',
    structure = SMILES('C#CC1OC([O])O1'),
    E0 = (27.2417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91479,0.0254097,-4.2907e-06,-3.71871e-09,9.97046e-13,3257.83,13.2757], Tmin=(100,'K'), Tmax=(2117.56,'K')), NASAPolynomial(coeffs=[19.2905,0.0141298,-8.88408e-06,1.68909e-09,-1.09116e-13,-7236.82,-81.8436], Tmin=(2117.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.2417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(OCOJ)"""),
)

species(
    label = 'C#CC1OC(O[O])O1(22525)',
    structure = SMILES('C#CC1OC(O[O])O1'),
    E0 = (23.3535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74542,0.0550015,-5.98803e-05,3.96496e-08,-1.15809e-11,2885.2,22.171], Tmin=(100,'K'), Tmax=(802.708,'K')), NASAPolynomial(coeffs=[6.45778,0.0315193,-1.59999e-05,3.20607e-09,-2.30725e-13,2128.67,0.472177], Tmin=(802.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.3535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(ROOJ)"""),
)

species(
    label = 'C#CC1OOC([O])O1(24740)',
    structure = SMILES('C#CC1OOC([O])O1'),
    E0 = (9.54002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65137,0.0463836,-3.54595e-05,1.56643e-08,-2.88598e-12,1236.37,22.7324], Tmin=(100,'K'), Tmax=(1285.53,'K')), NASAPolynomial(coeffs=[9.32417,0.0225089,-7.6013e-06,1.21704e-09,-7.63379e-14,-736.327,-16.2115], Tmin=(1285.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.54002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(124trioxolane) + radical(OCOJ)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[CH]OC([O])O[O](8400)',
    structure = SMILES('[CH]OC([O])O[O]'),
    E0 = (331.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,180,180,2106.73,2106.98,2108.16,2108.18],'cm^-1')),
        HinderedRotor(inertia=(0.129957,'amu*angstrom^2'), symmetry=1, barrier=(2.98798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131221,'amu*angstrom^2'), symmetry=1, barrier=(3.01702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129936,'amu*angstrom^2'), symmetry=1, barrier=(2.98749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05365,0.0655934,-0.000166056,1.93586e-07,-7.82481e-11,39971.1,21.8361], Tmin=(100,'K'), Tmax=(859.228,'K')), NASAPolynomial(coeffs=[-5.67417,0.0431564,-2.49126e-05,4.95398e-09,-3.43084e-13,43455.3,70.4934], Tmin=(859.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(OCOJ) + radical(CH2_triplet)"""),
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
    E0 = (216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (342.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (328.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (374.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (381.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (410.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (335.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (345.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (249.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (498.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (584.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (269.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (444.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (667.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (415.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (298.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (242.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (279.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (241.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (696.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (373.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (404.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (304.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (224.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (530.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (468.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (774.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (831.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (335.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (473.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (593.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (466.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (547.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (247.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (347.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (305.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (322.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (224.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (223.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (918.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[O]OC=O(5472)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=[C]C1OC(O[O])O1(24712)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(126.018,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 123.5 to 126.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=[C]C1OOC([O])O1(24713)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.51e+10,'s^-1'), n=0, Ea=(112.205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 107.8 to 112.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]=C=COC(=O)O[O](24714)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(48.8181,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 46.5 to 48.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC=O(5472)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.3e+11,'cm^3/(mol*s)'), n=0, Ea=(141.944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;O_rad/OneDe] for rate rule [CO-NdH_O;O_rad/OneDe]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 140.4 to 141.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', '[CH]=C=COC=O(18779)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.49e-08,'m^3/(mol*s)'), n=3.486, Ea=(224.094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-NdH_O;OJ] for rate rule [CO-NdH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 221.1 to 224.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=C=CO[C](O)O[O](24715)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH][C]=COC(=O)OO(24716)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C#COC([O])O[O](24717)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C=[C]OC(O)O[O](24718)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(222.678,'s^-1'), n=2.70078, Ea=(107.007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS_OCs;Y_rad_out;XH_out] for rate rule [R4H_SSS_OCs;Cd_rad_out_double;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C=[C]OC([O])OO(24719)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Cd_rad_out_double;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['C=[C][CH]OC(=O)O[O](24720)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][CH]O[O](8201)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC([O])[O](4172)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.43153e+08,'m^3/(mol*s)'), n=0.0716491, Ea=(15.4197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_allenic;O_rad/NonDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', '[CH][C]=COC=O(19586)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH][C]=COC(=O)O[O](24721)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]=C=[C]OC([O])O[O](24722)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[O]OC([O])OC1[C]=C1(24723)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(198.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 197.6 to 198.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]C1=COC(O[O])O1(24574)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.66666e+10,'s^-1'), n=0.302034, Ea=(82.5645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]C1=COC([O])OO1(24724)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(684195,'s^-1'), n=1.09371, Ea=(25.8736,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=C=COC(=O)OO(22526)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['C=C=COC(=O)O[O](24725)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]CO[C]([O])O[O](24726)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][CH][CH]OC(=O)O[O](24727)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['O2(S)(5486)', '[CH]=C=COC=O(18779)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['O(T)(63)', '[CH]=C=COC1OO1(24728)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(188.03,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['O(T)(63)', '[O]C1OC=C=CO1(24729)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.23936e+08,'s^-1'), n=0.739184, Ea=(87.8828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R5OO;Y_rad_intra;OO] for rate rule [R5OO;Cd_pri_rad_in;OOJ]
Euclidian distance = 2.2360679775
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=C=COC1OOO1(24730)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[CH]=C(C=O)C([O])O[O](22519)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(T)(63)', '[CH]=C=COC([O])[O](24731)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(T)(63)', '[CH]=C=CO[CH]O[O](19776)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/O2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(8)', '[C]=C=COC([O])O[O](24732)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O][CH]O[O](8201)', 'C#CC=O(21959)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.07e+06,'cm^3/(mol*s)'), n=2.43, Ea=(22.5936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CtH;YJ] for rate rule [Od_CO-CtH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#CCO[C]([O])O[O](24733)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[C]#CCOC([O])O[O](24734)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_1H] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]=C=COC(O)O[O](24735)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Ct_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]=C=COC([O])OO(24736)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(364667,'s^-1'), n=1.22214, Ea=(79.2357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;XH_out] for rate rule [R7HJ_2;Ct_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[O]OC1O[CH][C]=CO1(24608)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['[O]C1OC=[C][CH]OO1(24737)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['C#CCOC(=O)O[O](24738)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['O(T)(63)', 'C#CC1OC([O])O1(24739)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.58279e+10,'s^-1'), n=0.53, Ea=(106.301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;C_sec_rad_intra;OOJ] + [R3OO;C_sec_rad_intra;OO] for rate rule [R3OO;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['C#CC1OC(O[O])O1(22525)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C=COC([O])O[O](22517)'],
    products = ['C#CC1OOC([O])O1(24740)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R5_SSSS;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[C]#C(5143)', '[CH]OC([O])O[O](8400)'],
    products = ['[CH]=C=COC([O])O[O](22517)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

network(
    label = '4818',
    isomers = [
        '[CH]=C=COC([O])O[O](22517)',
    ],
    reactants = [
        ('[O]OC=O(5472)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4818',
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


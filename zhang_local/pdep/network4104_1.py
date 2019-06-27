species(
    label = '[CH]C(=C)C([CH2])CO[O](18660)',
    structure = SMILES('[CH]C(=C)C([CH2])CO[O]'),
    E0 = (556.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258411,0.0844173,-8.59158e-05,5.17606e-08,-1.31072e-11,67040.8,33.3871], Tmin=(100,'K'), Tmax=(943.563,'K')), NASAPolynomial(coeffs=[10.4818,0.0410773,-1.70164e-05,3.07959e-09,-2.08861e-13,65111.6,-15.3407], Tmin=(943.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
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
    label = 'C#C[CH2](17441)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,1131.03,1132.16,1135.9],'cm^-1')),
        HinderedRotor(inertia=(0.154206,'amu*angstrom^2'), symmetry=1, barrier=(3.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62061e-06,-1.82973e-08,7.68649e-12,39535.3,8.27851], Tmin=(100,'K'), Tmax=(960.555,'K')), NASAPolynomial(coeffs=[6.38511,0.00814486,-2.78734e-06,4.95348e-10,-3.50148e-14,38483.6,-8.79383], Tmin=(960.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = '[CH]C1([CH2])CC1CO[O](20981)',
    structure = SMILES('[CH]C1([CH2])CC1CO[O]'),
    E0 = (657.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242402,0.0714226,-4.63864e-05,6.06763e-09,3.43709e-12,79222.1,30.8238], Tmin=(100,'K'), Tmax=(1041.6,'K')), NASAPolynomial(coeffs=[17.6049,0.0254477,-9.98972e-06,1.85276e-09,-1.30859e-13,74482.2,-59.0384], Tmin=(1041.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Neopentyl) + radical(ROOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])OOCC1[CH2](20982)',
    structure = SMILES('[CH]C1([CH2])OOCC1[CH2]'),
    E0 = (631.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38633,0.0860854,-8.10576e-05,3.88809e-08,-6.87646e-12,76213.6,35.5753], Tmin=(100,'K'), Tmax=(1694.08,'K')), NASAPolynomial(coeffs=[18.4006,0.0174497,-8.80223e-07,-3.07429e-10,3.35868e-14,72654.3,-61.0342], Tmin=(1694.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(Isobutyl) + radical(CJCOOH) + radical(CCJ2_triplet)"""),
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
    label = '[CH]C(=C)C(=C)CO[O](20983)',
    structure = SMILES('[CH]C(=C)C(=C)CO[O]'),
    E0 = (459.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272572,0.0818734,-7.95608e-05,4.27494e-08,-9.46493e-12,55435.3,28.964], Tmin=(100,'K'), Tmax=(1079,'K')), NASAPolynomial(coeffs=[12.8515,0.0352416,-1.47341e-05,2.69574e-09,-1.84613e-13,52720.8,-32.6782], Tmin=(1079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
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
    label = '[CH]C(=C)C=C(17773)',
    structure = SMILES('[CH]C(=C)C=C'),
    E0 = (427.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.11235,'amu*angstrom^2'), symmetry=1, barrier=(48.5671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11944,'amu*angstrom^2'), symmetry=1, barrier=(48.7301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85297,0.0366833,7.83676e-06,-3.70295e-08,1.71666e-11,51444.9,17.4108], Tmin=(100,'K'), Tmax=(958.454,'K')), NASAPolynomial(coeffs=[11.8213,0.0206749,-7.16364e-06,1.2643e-09,-8.87592e-14,48358.5,-36.3899], Tmin=(958.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=C(18825)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.264,228.889,229.07],'cm^-1')),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27541,0.0127954,9.49515e-06,-1.56026e-08,5.42938e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.31,'K')), NASAPolynomial(coeffs=[4.18965,0.0168435,-6.77763e-06,1.22218e-09,-8.33556e-14,73336.3,4.89309], Tmin=(1063.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=C(C)CO[O](20984)',
    structure = SMILES('[CH]C([CH2])=C(C)CO[O]'),
    E0 = (484.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371927,0.0847337,-8.89715e-05,5.76244e-08,-1.61224e-11,58403.7,30.1888], Tmin=(100,'K'), Tmax=(847.4,'K')), NASAPolynomial(coeffs=[8.54535,0.0461515,-2.06746e-05,3.89249e-09,-2.70011e-13,57018.5,-7.88964], Tmin=(847.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(AllylJ2_triplet) + radical(Allyl_P) + radical(ROOJ)"""),
)

species(
    label = '[CH]C(=C)C([CH2])[CH]OO(20985)',
    structure = SMILES('[CH]C(=C)C([CH2])[CH]OO'),
    E0 = (592.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.049769,0.090084,-9.54546e-05,5.87012e-08,-1.50921e-11,71446.7,33.5301], Tmin=(100,'K'), Tmax=(931.011,'K')), NASAPolynomial(coeffs=[11.2157,0.0421108,-1.81627e-05,3.35514e-09,-2.30338e-13,69367.6,-19.5407], Tmin=(931.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(CCsJOOH)"""),
)

species(
    label = '[CH]C(=C)C(C)[CH]O[O](20986)',
    structure = SMILES('[CH]C(=C)C(C)[CH]O[O]'),
    E0 = (539.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386705,0.0836975,-8.36408e-05,5.01761e-08,-1.29402e-11,65049.9,31.301], Tmin=(100,'K'), Tmax=(915.713,'K')), NASAPolynomial(coeffs=[9.30358,0.0447468,-1.98366e-05,3.72443e-09,-2.58309e-13,63416.9,-10.9326], Tmin=(915.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(CCsJOOH) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=C([CH2])COO(20987)',
    structure = SMILES('[CH]C([CH2])=C([CH2])COO'),
    E0 = (484.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0980866,0.0861121,-7.73343e-05,3.76499e-08,-7.62973e-12,58356.8,30.6806], Tmin=(100,'K'), Tmax=(1160.83,'K')), NASAPolynomial(coeffs=[13.5381,0.0397998,-1.74899e-05,3.2808e-09,-2.27794e-13,55236.5,-36.1639], Tmin=(1160.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=[CH])C(C)CO[O](20988)',
    structure = SMILES('[CH]C(=[CH])C(C)CO[O]'),
    E0 = (598.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.302299,0.0846082,-8.41149e-05,4.86801e-08,-1.19288e-11,72091.4,31.1616], Tmin=(100,'K'), Tmax=(967.185,'K')), NASAPolynomial(coeffs=[10.5307,0.0423056,-1.85072e-05,3.4569e-09,-2.39223e-13,70112.9,-17.8432], Tmin=(967.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C([CH2])COO(20989)',
    structure = SMILES('[CH]C(=[CH])C([CH2])COO'),
    E0 = (651.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0457452,0.091118,-9.63274e-05,5.76854e-08,-1.42742e-11,78488.7,33.431], Tmin=(100,'K'), Tmax=(970.525,'K')), NASAPolynomial(coeffs=[12.4319,0.0396929,-1.68487e-05,3.09163e-09,-2.1161e-13,76066.7,-26.3933], Tmin=(970.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]C(=C)C([CH2])[CH2](17727)',
    structure = SMILES('[CH]C(=C)C([CH2])[CH2]'),
    E0 = (697.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,623.021,623.022,623.022,623.023],'cm^-1')),
        HinderedRotor(inertia=(0.200176,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200177,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200176,'amu*angstrom^2'), symmetry=1, barrier=(55.1376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200175,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.85,'J/mol'), sigma=(5.90911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.46 K, Pc=35.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21023,0.0555238,-2.63652e-05,-4.35895e-09,6.27756e-12,84020.2,26.8273], Tmin=(100,'K'), Tmax=(923.387,'K')), NASAPolynomial(coeffs=[10.335,0.0309114,-1.06121e-05,1.76031e-09,-1.15182e-13,81699.2,-19.9103], Tmin=(923.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(697.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C([CH2])=C[CH2](18836)',
    structure = SMILES('[CH]C([CH2])=C[CH2]'),
    E0 = (604.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,394.679,394.68,394.681],'cm^-1')),
        HinderedRotor(inertia=(0.458026,'amu*angstrom^2'), symmetry=1, barrier=(50.6299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458023,'amu*angstrom^2'), symmetry=1, barrier=(50.6299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458025,'amu*angstrom^2'), symmetry=1, barrier=(50.6299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81611,0.0388094,1.79234e-07,-2.5992e-08,1.22401e-11,72818.5,19.2374], Tmin=(100,'K'), Tmax=(992.533,'K')), NASAPolynomial(coeffs=[10.7608,0.0238542,-9.09622e-06,1.64932e-09,-1.15256e-13,70004,-29.0822], Tmin=(992.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(604.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])=C([CH2])CO[O](20990)',
    structure = SMILES('[CH]C([CH2])=C([CH2])CO[O]'),
    E0 = (636.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,492.5,1135,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372411,0.0823154,-8.10292e-05,4.54087e-08,-1.06924e-11,76627,30.4615], Tmin=(100,'K'), Tmax=(1007.18,'K')), NASAPolynomial(coeffs=[11.1441,0.0395354,-1.73156e-05,3.23515e-09,-2.23986e-13,74457.3,-21.5824], Tmin=(1007.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
)

species(
    label = '[CH]C(=C)C([CH2])[CH]O[O](20991)',
    structure = SMILES('[CH]C(=C)C([CH2])[CH]O[O]'),
    E0 = (744.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156732,0.0884539,-0.000107568,7.84662e-08,-2.37232e-11,89723.9,33.8976], Tmin=(100,'K'), Tmax=(821.232,'K')), NASAPolynomial(coeffs=[9.82831,0.040098,-1.69643e-05,3.06474e-09,-2.06049e-13,88177.4,-10.6011], Tmin=(821.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH]C(=[CH])C([CH2])CO[O](20992)',
    structure = SMILES('[CH]C(=[CH])C([CH2])CO[O]'),
    E0 = (803.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105127,0.0889072,-0.000106103,7.39581e-08,-2.11981e-11,96764,33.6452], Tmin=(100,'K'), Tmax=(846.363,'K')), NASAPolynomial(coeffs=[10.8818,0.0379723,-1.58256e-05,2.84381e-09,-1.90928e-13,94939.9,-16.5477], Tmin=(846.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1CCC1CO[O](20993)',
    structure = SMILES('[CH][C]1CCC1CO[O]'),
    E0 = (645.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14974,0.0628971,-4.28278e-05,1.52826e-08,-2.3152e-12,77743.6,29.6875], Tmin=(100,'K'), Tmax=(1462.87,'K')), NASAPolynomial(coeffs=[10.845,0.0363869,-1.56449e-05,2.89472e-09,-1.98161e-13,74907,-20.7747], Tmin=(1462.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCJ2_triplet) + radical(ROOJ) + radical(Tertalkyl)"""),
)

species(
    label = '[CH][C]1COOCC1[CH2](20994)',
    structure = SMILES('[CH][C]1COOCC1[CH2]'),
    E0 = (622.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37751,0.0384929,5.38769e-05,-1.07501e-07,4.86072e-11,74935,29.5853], Tmin=(100,'K'), Tmax=(889.79,'K')), NASAPolynomial(coeffs=[17.0002,0.0200618,-2.376e-06,7.29335e-11,2.23773e-16,70104.3,-55.4838], Tmin=(889.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(C2CJCOOH) + radical(CCJ2_triplet) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(=C)C(=C)COO(20995)',
    structure = SMILES('[CH]C(=C)C(=C)COO'),
    E0 = (307.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0664894,0.0863789,-7.81004e-05,3.75738e-08,-7.38338e-12,37168,29.4191], Tmin=(100,'K'), Tmax=(1210.88,'K')), NASAPolynomial(coeffs=[15.406,0.0352674,-1.47852e-05,2.71483e-09,-1.86378e-13,33420.9,-48.1875], Tmin=(1210.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[C]([CH2])CO[O](20996)',
    structure = SMILES('[CH]C([CH2])[C]([CH2])CO[O]'),
    E0 = (909.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,199.681,1048.08,1296.61,1887.36,2135.25],'cm^-1')),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124127,'amu*angstrom^2'), symmetry=1, barrier=(2.89416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666872,0.0765954,-7.98287e-05,4.52373e-08,-8.80053e-12,109531,39.1543], Tmin=(100,'K'), Tmax=(725.487,'K')), NASAPolynomial(coeffs=[9.81883,0.0354532,-1.40286e-05,2.47469e-09,-1.64993e-13,107958,-3.75152], Tmin=(725.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(909.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Isobutyl) + radical(ROOJ) + radical(C2CJCOOH)"""),
)

species(
    label = '[CH]C([CH2])C([CH2])[CH]O[O](20997)',
    structure = SMILES('[CH]C([CH2])C([CH2])[CH]O[O]'),
    E0 = (907.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.136351,0.0955558,-0.000131393,1.01163e-07,-3.07874e-11,109279,37.3218], Tmin=(100,'K'), Tmax=(903.165,'K')), NASAPolynomial(coeffs=[11.5194,0.0342686,-1.35534e-05,2.33166e-09,-1.50579e-13,107567,-15.5412], Tmin=(903.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(ROOJ) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCsJOOH)"""),
)

species(
    label = '[CH][C](C)C([CH2])[CH]O[O](20998)',
    structure = SMILES('[CH][C](C)C([CH2])[CH]O[O]'),
    E0 = (887.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,218.267,812.164,1082.89,1394.22,1697.32],'cm^-1')),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137704,'amu*angstrom^2'), symmetry=1, barrier=(3.56292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0213244,0.100769,-0.00016417,1.49097e-07,-5.18422e-11,106903,36.2609], Tmin=(100,'K'), Tmax=(880.934,'K')), NASAPolynomial(coeffs=[6.17871,0.0443798,-2.00729e-05,3.66235e-09,-2.43427e-13,106906,13.354], Tmin=(880.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(ROOJ) + radical(Isobutyl) + radical(CCsJOOH) + radical(Tertalkyl)"""),
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
    label = '[CH]C(=C)C1COC1(20999)',
    structure = SMILES('[CH]C(=C)C1COC1'),
    E0 = (307.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44057,0.0426723,3.17824e-05,-7.45898e-08,3.41022e-11,37095.9,22.5527], Tmin=(100,'K'), Tmax=(891.904,'K')), NASAPolynomial(coeffs=[12.6229,0.0293031,-7.5916e-06,1.07762e-09,-6.74748e-14,33638.2,-38.3174], Tmin=(891.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Oxetane) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=CCCO[O](21000)',
    structure = SMILES('[CH]C([CH2])=CCCO[O]'),
    E0 = (494.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435499,0.0813536,-7.43556e-05,3.92059e-08,-8.85149e-12,59572.1,31.0586], Tmin=(100,'K'), Tmax=(1035.65,'K')), NASAPolynomial(coeffs=[10.3335,0.0431242,-1.89851e-05,3.56263e-09,-2.4735e-13,57521.9,-17.0401], Tmin=(1035.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C[CH]CO[O](18713)',
    structure = SMILES('[CH]C(=C)C[CH]CO[O]'),
    E0 = (555.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294144,0.0870235,-9.88667e-05,7.0375e-08,-2.15341e-11,66930.5,34.032], Tmin=(100,'K'), Tmax=(780.724,'K')), NASAPolynomial(coeffs=[8.19305,0.0465552,-2.11181e-05,3.98708e-09,-2.76409e-13,65697.1,-2.1207], Tmin=(780.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCJCOOH) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
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
    label = '[CH]C(=C)C([CH2])=C(18076)',
    structure = SMILES('[CH]C(=C)C([CH2])=C'),
    E0 = (540.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,195.967,195.967,195.967,195.967],'cm^-1')),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3249.64,'J/mol'), sigma=(5.69809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.59 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01928,0.0540671,-1.36977e-05,-2.39883e-08,1.40947e-11,65175.7,20.5913], Tmin=(100,'K'), Tmax=(963.045,'K')), NASAPolynomial(coeffs=[14.94,0.0233813,-8.16493e-06,1.43762e-09,-1.00384e-13,61236.2,-52.5766], Tmin=(963.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C1COOC1(21001)',
    structure = SMILES('[CH]C(=C)C1COOC1'),
    E0 = (279.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854317,0.0524202,2.6618e-05,-7.8968e-08,3.74596e-11,33740.1,26.1356], Tmin=(100,'K'), Tmax=(900.644,'K')), NASAPolynomial(coeffs=[16.5581,0.0275202,-6.59955e-06,9.04872e-10,-5.73991e-14,29092.6,-58.0803], Tmin=(900.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(12dioxolane) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C([CH2])C[O](18167)',
    structure = SMILES('[CH]C(=C)C([CH2])C[O]'),
    E0 = (558.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08553,0.0676325,-5.95345e-05,3.49777e-08,-9.43518e-12,67273.8,29.763], Tmin=(100,'K'), Tmax=(854.509,'K')), NASAPolynomial(coeffs=[5.98482,0.044698,-1.92743e-05,3.56679e-09,-2.45161e-13,66436.5,6.89717], Tmin=(854.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(CCOJ) + radical(Isobutyl)"""),
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
    label = '[CH]C(=C)[CH]CO[O](21002)',
    structure = SMILES('[CH]C(=C)[CH]CO[O]'),
    E0 = (495.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.818968,0.0727533,-7.19262e-05,4.03807e-08,-9.55505e-12,59732,24.7951], Tmin=(100,'K'), Tmax=(1000.14,'K')), NASAPolynomial(coeffs=[10.1717,0.035348,-1.58266e-05,2.98658e-09,-2.07924e-13,57861.1,-20.3278], Tmin=(1000.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(AllylJ2_triplet) + radical(ROOJ)"""),
)

species(
    label = '[CH]C(=C)C([CH])CO[O](21003)',
    structure = SMILES('[CH]C(=C)C([CH])CO[O]'),
    E0 = (799.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.300305,0.084425,-8.80083e-05,5.23807e-08,-1.30188e-11,96280.7,32.4622], Tmin=(100,'K'), Tmax=(960.761,'K')), NASAPolynomial(coeffs=[11.2033,0.039031,-1.7135e-05,3.20102e-09,-2.21451e-13,94185.7,-19.7016], Tmin=(960.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(799.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(AllylJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C#CC([CH2])CO[O](7729)',
    structure = SMILES('C#CC([CH2])CO[O]'),
    E0 = (385.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(3.62051,'amu*angstrom^2'), symmetry=1, barrier=(83.2427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141226,'amu*angstrom^2'), symmetry=1, barrier=(3.24707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141807,'amu*angstrom^2'), symmetry=1, barrier=(3.26043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.62318,'amu*angstrom^2'), symmetry=1, barrier=(83.304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84165,0.0711844,-9.19474e-05,6.66223e-08,-1.90454e-11,46469.9,27.2268], Tmin=(100,'K'), Tmax=(952.029,'K')), NASAPolynomial(coeffs=[10.3765,0.0247813,-8.84299e-06,1.43036e-09,-8.86905e-14,44941.8,-16.795], Tmin=(952.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = '[CH]C(C)=C([CH2])CO[O](21004)',
    structure = SMILES('[CH]C(C)=C([CH2])CO[O]'),
    E0 = (484.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371927,0.0847337,-8.89715e-05,5.76244e-08,-1.61224e-11,58403.7,30.1888], Tmin=(100,'K'), Tmax=(847.4,'K')), NASAPolynomial(coeffs=[8.54535,0.0461515,-2.06746e-05,3.89249e-09,-2.70011e-13,57018.5,-7.88964], Tmin=(847.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])=C([CH2])CO[O](18658)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])CO[O]'),
    E0 = (416.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,492.5,1135,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2425.72],'cm^-1')),
        HinderedRotor(inertia=(0.0322788,'amu*angstrom^2'), symmetry=1, barrier=(1.88062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46125,'amu*angstrom^2'), symmetry=1, barrier=(81.5766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575844,'amu*angstrom^2'), symmetry=1, barrier=(33.3605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362448,'amu*angstrom^2'), symmetry=1, barrier=(20.253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40082,'amu*angstrom^2'), symmetry=1, barrier=(81.5821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242264,0.0815256,-7.89534e-05,4.09729e-08,-8.65426e-12,50272.9,30.093], Tmin=(100,'K'), Tmax=(1133.95,'K')), NASAPolynomial(coeffs=[14.2644,0.0320622,-1.35221e-05,2.50439e-09,-1.73095e-13,47092.9,-39.3181], Tmin=(1133.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C)C([CH2])[CH]O[O](21005)',
    structure = SMILES('[CH]=C(C)C([CH2])[CH]O[O]'),
    E0 = (621.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.178321,'amu*angstrom^2'), symmetry=1, barrier=(4.09996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177008,'amu*angstrom^2'), symmetry=1, barrier=(4.06975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17742,'amu*angstrom^2'), symmetry=1, barrier=(4.07923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180849,'amu*angstrom^2'), symmetry=1, barrier=(4.15808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.80968,'amu*angstrom^2'), symmetry=1, barrier=(64.6001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0283522,0.0927463,-0.000127295,1.00319e-07,-3.17598e-11,74862.8,33.6498], Tmin=(100,'K'), Tmax=(853.664,'K')), NASAPolynomial(coeffs=[10.3775,0.0367366,-1.56707e-05,2.83125e-09,-1.89359e-13,73369.7,-13.0373], Tmin=(853.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ) + radical(Isobutyl) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])[CH]O[O](18659)',
    structure = SMILES('[CH2]C(=C)C([CH2])[CH]O[O]'),
    E0 = (525.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,314.35,314.416],'cm^-1')),
        HinderedRotor(inertia=(0.00170614,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0941533,'amu*angstrom^2'), symmetry=1, barrier=(6.60308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06105,'amu*angstrom^2'), symmetry=1, barrier=(74.3724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012252,'amu*angstrom^2'), symmetry=1, barrier=(74.373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06051,'amu*angstrom^2'), symmetry=1, barrier=(74.3722,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211265,0.0853466,-9.67595e-05,6.18451e-08,-1.61103e-11,63361.9,32.1833], Tmin=(100,'K'), Tmax=(930.42,'K')), NASAPolynomial(coeffs=[12.1564,0.0339929,-1.39688e-05,2.52405e-09,-1.71023e-13,61139.1,-24.5838], Tmin=(930.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCsJOOH) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(C)C(=C)CO[O](21006)',
    structure = SMILES('[CH]=C(C)C(=C)CO[O]'),
    E0 = (336.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229797,0.0850468,-9.48887e-05,5.8264e-08,-1.45432e-11,40570.7,28.4167], Tmin=(100,'K'), Tmax=(968.951,'K')), NASAPolynomial(coeffs=[12.8916,0.0327764,-1.39706e-05,2.58983e-09,-1.78655e-13,38117,-32.27], Tmin=(968.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=C)C(=C)CO[O](18648)',
    structure = SMILES('[CH2]C(=C)C(=C)CO[O]'),
    E0 = (240.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,263.281,263.997],'cm^-1')),
        HinderedRotor(inertia=(0.226491,'amu*angstrom^2'), symmetry=1, barrier=(11.1505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225393,'amu*angstrom^2'), symmetry=1, barrier=(11.148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43006,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41176,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0817632,0.0817437,-7.95512e-05,4.06895e-08,-8.32662e-12,29083.9,28.1239], Tmin=(100,'K'), Tmax=(1181.93,'K')), NASAPolynomial(coeffs=[16.0718,0.0276283,-1.08726e-05,1.95115e-09,-1.32717e-13,25304.1,-51.6917], Tmin=(1181.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C1COCC1[CH2](21007)',
    structure = SMILES('[CH]=C1COCC1[CH2]'),
    E0 = (318.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82337,0.0316781,4.59511e-05,-8.22026e-08,3.4124e-11,38445.9,25.0249], Tmin=(100,'K'), Tmax=(944.149,'K')), NASAPolynomial(coeffs=[13.1071,0.0246233,-7.58143e-06,1.31057e-09,-9.38224e-14,34498.9,-38.3824], Tmin=(944.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1=COCC1[CH2](21008)',
    structure = SMILES('[CH2]C1=COCC1[CH2]'),
    E0 = (157.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37414,0.0333073,7.07204e-05,-1.33432e-07,6.01432e-11,19030.5,19.8937], Tmin=(100,'K'), Tmax=(893.229,'K')), NASAPolynomial(coeffs=[21.6967,0.00786555,3.34079e-06,-9.66185e-10,6.83211e-14,12784.4,-90.4973], Tmin=(893.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]CC([CH2])CO[O](18633)',
    structure = SMILES('[CH]=[C]CC([CH2])CO[O]'),
    E0 = (680.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.149134,'amu*angstrom^2'), symmetry=1, barrier=(3.42889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149193,'amu*angstrom^2'), symmetry=1, barrier=(3.43025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149337,'amu*angstrom^2'), symmetry=1, barrier=(3.43356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149022,'amu*angstrom^2'), symmetry=1, barrier=(3.4263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755,'amu*angstrom^2'), symmetry=1, barrier=(17.3589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3940.69,'J/mol'), sigma=(6.78976,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.53 K, Pc=28.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.128683,0.0903961,-0.000122702,9.64656e-08,-3.06141e-11,81989.1,35.0673], Tmin=(100,'K'), Tmax=(845.331,'K')), NASAPolynomial(coeffs=[10.0688,0.0367428,-1.57534e-05,2.85988e-09,-1.92046e-13,80545,-9.81933], Tmin=(845.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CO[O])C1=CC1(21009)',
    structure = SMILES('[CH2]C(CO[O])C1=CC1'),
    E0 = (449.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.304176,0.0764841,-7.15009e-05,3.61216e-08,-7.35955e-12,54255.2,31.6939], Tmin=(100,'K'), Tmax=(1184.17,'K')), NASAPolynomial(coeffs=[14.5336,0.0284182,-1.06149e-05,1.84365e-09,-1.22807e-13,50885.2,-39.3604], Tmin=(1184.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1CCC1CO[O](20977)',
    structure = SMILES('[CH]=C1CCC1CO[O]'),
    E0 = (379.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.2691,-0.00619483,0.00011626,-1.14696e-07,2.61338e-11,45232.5,-17.0869], Tmin=(100,'K'), Tmax=(1762.56,'K')), NASAPolynomial(coeffs=[83.721,0.0204118,-6.89573e-05,1.6851e-08,-1.24662e-12,-9980.71,-491.865], Tmin=(1762.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C1[CH]CC1CO[O](21010)',
    structure = SMILES('C=C1[CH]CC1CO[O]'),
    E0 = (273.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.0779,-0.00989208,0.000125928,-1.2201e-07,2.79037e-11,32499.5,-17.8807], Tmin=(100,'K'), Tmax=(1739.43,'K')), NASAPolynomial(coeffs=[79.1796,0.0262705,-7.14923e-05,1.73681e-08,-1.28556e-12,-20354.3,-467.957], Tmin=(1739.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(ROOJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C1COOCC1[CH2](21011)',
    structure = SMILES('[CH]=C1COOCC1[CH2]'),
    E0 = (348.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903332,0.0500664,1.93596e-05,-5.81018e-08,2.41567e-11,41998.3,25.5993], Tmin=(100,'K'), Tmax=(1030.76,'K')), NASAPolynomial(coeffs=[16.8469,0.0303327,-1.3243e-05,2.6447e-09,-1.95867e-13,36473,-62.6608], Tmin=(1030.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1COO[CH]C1=C(21012)',
    structure = SMILES('[CH2]C1COO[CH]C1=C'),
    E0 = (218.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14275,0.0360296,7.38739e-05,-1.22924e-07,4.88497e-11,26386.4,24.9314], Tmin=(100,'K'), Tmax=(985.995,'K')), NASAPolynomial(coeffs=[20.2977,0.0247302,-9.96448e-06,2.07144e-09,-1.62675e-13,19381,-83.5799], Tmin=(985.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(Isobutyl) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=[C]C([CH2])CO[O](7740)',
    structure = SMILES('[CH]=[C]C([CH2])CO[O]'),
    E0 = (705.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.174848,'amu*angstrom^2'), symmetry=1, barrier=(4.0201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174878,'amu*angstrom^2'), symmetry=1, barrier=(4.02078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174937,'amu*angstrom^2'), symmetry=1, barrier=(4.02214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.719515,'amu*angstrom^2'), symmetry=1, barrier=(16.5431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79242,0.075123,-0.000106546,8.49227e-08,-2.68054e-11,84954.6,30.4739], Tmin=(100,'K'), Tmax=(877.865,'K')), NASAPolynomial(coeffs=[9.3267,0.0279161,-1.16676e-05,2.07354e-09,-1.36773e-13,83776.8,-7.76128], Tmin=(877.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]C(=C)C([CH2])CO[O](21013)',
    structure = SMILES('[C]C(=C)C([CH2])CO[O]'),
    E0 = (855.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,263.191,263.192,263.194],'cm^-1')),
        HinderedRotor(inertia=(0.15977,'amu*angstrom^2'), symmetry=1, barrier=(7.85383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243357,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159779,'amu*angstrom^2'), symmetry=1, barrier=(7.85389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24794,'amu*angstrom^2'), symmetry=1, barrier=(110.5,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.131692,0.0875115,-0.000109714,7.51101e-08,-2.05796e-11,102982,31.1841], Tmin=(100,'K'), Tmax=(892.402,'K')), NASAPolynomial(coeffs=[13.0931,0.0294168,-1.20686e-05,2.16665e-09,-1.45721e-13,100668,-29.8721], Tmin=(892.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(CJ3) + radical(Isobutyl) + radical(ROOJ)"""),
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
    E0 = (556.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (657.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (631.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (671.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (663.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (704.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (687.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (722.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (673.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (686.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (642.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (713.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (700.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (804.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (963.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (847.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (956.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1015.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (682.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (622.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (581.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (932.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (970.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (896.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (639.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (716.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (716.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (695.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (563.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (806.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (911.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1011.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (782.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (701.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (758.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (748.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (717.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (600.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (619.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (619.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (618.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (611.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (926.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (564.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (564.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (564.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (565.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (564.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (1121.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (1066.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['C=CCO[O](6082)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C1([CH2])CC1CO[O](20981)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+10,'s^-1'), n=0.35, Ea=(101.18,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 98.7 to 101.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C1([CH2])OOCC1[CH2](20982)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.51e+10,'s^-1'), n=0, Ea=(75.5267,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra_2H;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 73.0 to 75.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]C(=C)C(=C)CO[O](20983)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.79403,'m^3/(mol*s)'), n=1.96942, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]O[O](61)', '[CH]C(=C)C=C(17773)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00203367,'m^3/(mol*s)'), n=2.41, Ea=(36.063,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;CsJ] for rate rule [Cds-CdH_Cds-HH;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CCO[O](6082)', '[CH][C]=C(18825)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.540096,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C([CH2])=C(C)CO[O](20984)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(=C)C([CH2])[CH]OO(20985)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.66e+08,'s^-1'), n=1.28, Ea=(166.272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 241 used for R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(=C)C(C)[CH]O[O](20986)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C([CH2])=C([CH2])COO(20987)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17136e+07,'s^-1'), n=1.54267, Ea=(130.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS_OCs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C(=[CH])C(C)CO[O](20988)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(=[CH])C([CH2])COO(20989)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(26171.5,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O2(2)', '[CH]C(=C)C([CH2])[CH2](17727)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.66237e+08,'m^3/(mol*s)'), n=-0.783071, Ea=(11.596,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for O2_birad;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]O[O](61)', '[CH]C([CH2])=C[CH2](18836)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.88087e+06,'m^3/(mol*s)'), n=0.114385, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_pri_rad] for rate rule [Y_rad;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CO[O](6081)', '[CH][C]=C(18825)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]C([CH2])=C([CH2])CO[O](20990)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]C(=C)C([CH2])[CH]O[O](20991)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]C(=[CH])C([CH2])CO[O](20992)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.15742e+08,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH][C]1CCC1CO[O](20993)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH][C]1COOCC1[CH2](20994)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(684195,'s^-1'), n=1.09371, Ea=(65.8052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 59.3 to 65.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(=C)C(=C)COO(20995)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C([CH2])[C]([CH2])CO[O](20996)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C([CH2])C([CH2])[CH]O[O](20997)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][C](C)C([CH2])[CH]O[O](20998)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['O(T)(63)', '[CH]C(=C)C1COC1(20999)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C([CH2])=CCCO[O](21000)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(=C)C[CH]CO[O](18713)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[O]O(16)', '[CH]C(=C)C([CH2])=C(18076)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(=C)C1COOC1(21001)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(T)(63)', '[CH]C(=C)C([CH2])C[O](18167)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(T)(28)', '[CH]C(=C)[CH]CO[O](21002)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(8)', '[CH]C(=C)C([CH])CO[O](21003)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(T)(28)', 'C#CC([CH2])CO[O](7729)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CO[O](6081)', 'C#C[CH2](17441)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.0031273,'m^3/(mol*s)'), n=2.54618, Ea=(24.6367,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]C(C)=C([CH2])CO[O](21004)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH2]C([CH2])=C([CH2])CO[O](18658)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C(C)C([CH2])[CH]O[O](21005)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(83700,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 321 used for R4H_SSS;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH2]C(=C)C([CH2])[CH]O[O](18659)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]=C(C)C(=C)CO[O](21006)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH2]C(=C)C(=C)CO[O](18648)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['O(T)(63)', '[CH]=C1COCC1[CH2](21007)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;C_pri_rad_intra;OO] for rate rule [R4OO_SSS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['O(T)(63)', '[CH2]C1=COCC1[CH2](21008)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_DSS;Cd_pri_rad_in;OOJ]
Euclidian distance = 2.44948974278
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]CC([CH2])CO[O](18633)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH2]C(CO[O])C1=CC1(21009)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]=C1CCC1CO[O](20977)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['C=C1[CH]CC1CO[O](21010)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH]=C1COOCC1[CH2](21011)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R6_SSSSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSSSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    products = ['[CH2]C1COO[CH]C1=C(21012)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CH2(T)(28)', '[CH]=[C]C([CH2])CO[O](7740)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['H(8)', '[C]C(=C)C([CH2])CO[O](21013)'],
    products = ['[CH]C(=C)C([CH2])CO[O](18660)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '4104',
    isomers = [
        '[CH]C(=C)C([CH2])CO[O](18660)',
    ],
    reactants = [
        ('C=CCO[O](6082)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4104',
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


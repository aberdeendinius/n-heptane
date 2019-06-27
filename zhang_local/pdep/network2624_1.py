species(
    label = '[CH2]C([CH2])([O])C[CH]CC(9519)',
    structure = SMILES('[CH2]C([CH2])([O])C[CH]CC'),
    E0 = (444.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,382.871,780.132,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155166,'amu*angstrom^2'), symmetry=1, barrier=(3.56757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.697188,0.112076,-0.000156128,1.30476e-07,-4.43231e-11,53640.4,34.9339], Tmin=(100,'K'), Tmax=(821.085,'K')), NASAPolynomial(coeffs=[9.0434,0.0521195,-2.37544e-05,4.45048e-09,-3.04597e-13,52462.3,-7.57203], Tmin=(821.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=CCC(36)',
    structure = SMILES('C=CCC'),
    E0 = (-16.5218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,385.42],'cm^-1')),
        HinderedRotor(inertia=(0.10941,'amu*angstrom^2'), symmetry=1, barrier=(11.0707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106381,'amu*angstrom^2'), symmetry=1, barrier=(11.1092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64174,0.0203702,3.0581e-05,-4.85559e-08,1.85232e-11,-1929.81,14.9537], Tmin=(100,'K'), Tmax=(991.794,'K')), NASAPolynomial(coeffs=[7.80835,0.0229241,-8.65887e-06,1.6005e-09,-1.13877e-13,-4105.1,-15.7294], Tmin=(991.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.5218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C([CH2])([O])C=CCC(10454)',
    structure = SMILES('[CH2]C([CH2])([O])C=CCC'),
    E0 = (359.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1572.82,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148028,'amu*angstrom^2'), symmetry=1, barrier=(3.40345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148028,'amu*angstrom^2'), symmetry=1, barrier=(3.40345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148028,'amu*angstrom^2'), symmetry=1, barrier=(3.40345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148028,'amu*angstrom^2'), symmetry=1, barrier=(3.40345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148028,'amu*angstrom^2'), symmetry=1, barrier=(3.40345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288139,0.0821158,-6.76601e-05,3.06869e-08,-5.87433e-12,43331.3,34.187], Tmin=(100,'K'), Tmax=(1210.73,'K')), NASAPolynomial(coeffs=[12.4093,0.0420698,-1.8046e-05,3.36767e-09,-2.33232e-13,40396.2,-26.6082], Tmin=(1210.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([CH2])([O])CC=CC(10455)',
    structure = SMILES('[CH2]C([CH2])([O])CC=CC'),
    E0 = (362.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1620.7,2888.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150258,'amu*angstrom^2'), symmetry=1, barrier=(3.45473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150258,'amu*angstrom^2'), symmetry=1, barrier=(3.45473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150258,'amu*angstrom^2'), symmetry=1, barrier=(3.45473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150258,'amu*angstrom^2'), symmetry=1, barrier=(3.45473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150258,'amu*angstrom^2'), symmetry=1, barrier=(3.45473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.35334,0.100934,-0.000115174,7.46576e-08,-2.00082e-11,43737.3,31.1873], Tmin=(100,'K'), Tmax=(897.958,'K')), NASAPolynomial(coeffs=[12.4874,0.0437351,-1.96276e-05,3.7232e-09,-2.59808e-13,41431.2,-29.3798], Tmin=(897.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C([O])C[CH]CC(950)',
    structure = SMILES('C=C([O])C[CH]CC'),
    E0 = (53.6223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,2699.61,2699.67],'cm^-1')),
        HinderedRotor(inertia=(0.206614,'amu*angstrom^2'), symmetry=1, barrier=(4.75047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231681,'amu*angstrom^2'), symmetry=1, barrier=(13.984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608166,'amu*angstrom^2'), symmetry=1, barrier=(13.9829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231711,'amu*angstrom^2'), symmetry=1, barrier=(13.9843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1185,0.0634505,-4.53405e-05,1.82342e-08,-3.18071e-12,6552.72,29.0296], Tmin=(100,'K'), Tmax=(1288.34,'K')), NASAPolynomial(coeffs=[9.34835,0.0378989,-1.55913e-05,2.84021e-09,-1.9357e-13,4432.14,-12.7598], Tmin=(1288.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.6223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC)"""),
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
    label = '[CH2]C(=C)C[CH]CC(6234)',
    structure = SMILES('[CH2]C(=C)C[CH]CC'),
    E0 = (242.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839832,0.0669513,-3.9312e-05,1.16119e-08,-1.42369e-12,29319.6,28.8165], Tmin=(100,'K'), Tmax=(1779.23,'K')), NASAPolynomial(coeffs=[13.2578,0.0390338,-1.5776e-05,2.79316e-09,-1.84568e-13,24900.7,-38.2481], Tmin=(1779.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CC(39)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1553.58],'cm^-1')),
        HinderedRotor(inertia=(0.00260974,'amu*angstrom^2'), symmetry=1, barrier=(4.49727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19501,'amu*angstrom^2'), symmetry=1, barrier=(4.48366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194146,'amu*angstrom^2'), symmetry=1, barrier=(4.4638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59221,0.0286125,-6.18652e-06,-3.09336e-09,1.2171e-12,30768.7,19.1946], Tmin=(100,'K'), Tmax=(1492.38,'K')), NASAPolynomial(coeffs=[6.22741,0.0257398,-1.0205e-05,1.78668e-09,-1.17174e-13,28918.5,-2.36199], Tmin=(1492.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH2])[O](10271)',
    structure = SMILES('[CH2][C]([CH2])[O]'),
    E0 = (537.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,278.503],'cm^-1')),
        HinderedRotor(inertia=(0.00215299,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939305,'amu*angstrom^2'), symmetry=1, barrier=(5.13965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0694,0.0447257,-6.5608e-05,5.12452e-08,-1.57124e-11,64674.3,16.4544], Tmin=(100,'K'), Tmax=(870.707,'K')), NASAPolynomial(coeffs=[8.27065,0.0130302,-5.47972e-06,9.76923e-10,-6.45299e-14,63716,-11.9064], Tmin=(870.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CJCO) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH3](11)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1408.13,1408.49,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48742e-06,-3.32748e-09,8.49957e-13,16285.6,0.351741], Tmin=(100,'K'), Tmax=(1337.63,'K')), NASAPolynomial(coeffs=[3.54146,0.00476787,-1.82148e-06,3.28877e-10,-2.22546e-14,16224,1.66035], Tmin=(1337.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH2]C([CH2])([O])CC=C(10187)',
    structure = SMILES('[CH2]C([CH2])([O])CC=C'),
    E0 = (398.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1570.56,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.145819,'amu*angstrom^2'), symmetry=1, barrier=(3.35267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145819,'amu*angstrom^2'), symmetry=1, barrier=(3.35267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145819,'amu*angstrom^2'), symmetry=1, barrier=(3.35267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145819,'amu*angstrom^2'), symmetry=1, barrier=(3.35267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216741,0.0863366,-9.83632e-05,6.14974e-08,-1.56531e-11,48051.5,27.3295], Tmin=(100,'K'), Tmax=(949.608,'K')), NASAPolynomial(coeffs=[12.7878,0.0333835,-1.47179e-05,2.77435e-09,-1.93175e-13,45664,-32.6686], Tmin=(949.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([CH2])([O])[CH]CCC(10423)',
    structure = SMILES('[CH2]C([CH2])([O])[CH]CCC'),
    E0 = (450.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.518901,0.103385,-0.000112008,6.80682e-08,-1.70629e-11,54293.7,34.2398], Tmin=(100,'K'), Tmax=(957.881,'K')), NASAPolynomial(coeffs=[13.6305,0.0442983,-1.9481e-05,3.67092e-09,-2.55638e-13,51583.1,-33.4139], Tmin=(957.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CC(C)2OJ) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C([CH2])([O])CC[CH]C(10456)',
    structure = SMILES('[CH2]C([CH2])([O])CC[CH]C'),
    E0 = (444.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723231,0.111453,-0.000150347,1.2121e-07,-4.00485e-11,53641,35.2002], Tmin=(100,'K'), Tmax=(810.428,'K')), NASAPolynomial(coeffs=[10.0955,0.0502408,-2.2588e-05,4.21684e-09,-2.8855e-13,52144.1,-13.1368], Tmin=(810.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([CH2])(O)[CH][CH]CC(10457)',
    structure = SMILES('[CH2]C([CH2])(O)[CH][CH]CC'),
    E0 = (415.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.728826,0.1109,-0.000149695,1.18647e-07,-3.83055e-11,50153.8,37.4843], Tmin=(100,'K'), Tmax=(820.6,'K')), NASAPolynomial(coeffs=[10.9593,0.0471472,-2.07668e-05,3.83703e-09,-2.60882e-13,48463.8,-15.2024], Tmin=(820.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCJCO) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(C)([O])[CH][CH]CC(10458)',
    structure = SMILES('[CH2]C(C)([O])[CH][CH]CC'),
    E0 = (432.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00278736,0.0937336,-0.000100601,6.72842e-08,-1.94067e-11,52161.4,35.7888], Tmin=(100,'K'), Tmax=(824.177,'K')), NASAPolynomial(coeffs=[8.81748,0.0509508,-2.27324e-05,4.29367e-09,-2.98576e-13,50708.5,-5.032], Tmin=(824.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(CCJCO) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCCC([CH2])([CH2])[O](10459)',
    structure = SMILES('[CH2]CCCC([CH2])([CH2])[O]'),
    E0 = (455.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,180,180,373.656,775.705,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154607,'amu*angstrom^2'), symmetry=1, barrier=(3.55471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.712148,0.109717,-0.000132274,9.11992e-08,-2.5909e-11,54941.4,34.2259], Tmin=(100,'K'), Tmax=(850.831,'K')), NASAPolynomial(coeffs=[12.7221,0.0465599,-2.09295e-05,3.95699e-09,-2.74966e-13,52655.3,-28.4166], Tmin=(850.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CC(C)2OJ) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C([CH2])(O)C[CH][CH]C(10460)',
    structure = SMILES('[CH2]C([CH2])(O)C[CH][CH]C'),
    E0 = (410.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.740117,0.116553,-0.000178825,1.58463e-07,-5.48317e-11,49492.8,37.7606], Tmin=(100,'K'), Tmax=(866.058,'K')), NASAPolynomial(coeffs=[7.28597,0.0533251,-2.40091e-05,4.41472e-09,-2.96415e-13,49083.6,5.85723], Tmin=(866.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(C)([O])C[CH][CH]C(10461)',
    structure = SMILES('[CH2]C(C)([O])C[CH][CH]C'),
    E0 = (427.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.191813,0.101776,-0.000139252,1.21341e-07,-4.29721e-11,51508.2,36.7085], Tmin=(100,'K'), Tmax=(833.689,'K')), NASAPolynomial(coeffs=[5.55319,0.0563906,-2.55309e-05,4.76333e-09,-3.24946e-13,51169.6,13.7514], Tmin=(833.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(RCCJC) + radical(CJC(C)2O) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])([CH2])O(10462)',
    structure = SMILES('[CH2]C[CH]CC([CH2])([CH2])O'),
    E0 = (420.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.908634,0.117143,-0.000169955,1.42146e-07,-4.74399e-11,50800.8,37.4177], Tmin=(100,'K'), Tmax=(844.906,'K')), NASAPolynomial(coeffs=[10.2708,0.0490008,-2.19651e-05,4.06126e-09,-2.74905e-13,49454.8,-11.4186], Tmin=(844.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(RCCJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])(C)[O](10463)',
    structure = SMILES('[CH2]C[CH]CC([CH2])(C)[O]'),
    E0 = (437.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.308744,0.101707,-0.000127843,1.01413e-07,-3.39129e-11,52814,36.1839], Tmin=(100,'K'), Tmax=(781.702,'K')), NASAPolynomial(coeffs=[8.36016,0.0523852,-2.36778e-05,4.4562e-09,-3.07358e-13,51610.3,-2.5338], Tmin=(781.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(RCCJ) + radical(CJC(C)2O) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C([CH2])([O])[CH][CH]CC(10464)',
    structure = SMILES('[CH2]C([CH2])([O])[CH][CH]CC'),
    E0 = (644.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,331.988,830.571,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153894,'amu*angstrom^2'), symmetry=1, barrier=(3.53832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.511506,0.106956,-0.000147767,1.21202e-07,-4.04939e-11,77677.4,36.7645], Tmin=(100,'K'), Tmax=(816.095,'K')), NASAPolynomial(coeffs=[9.64552,0.0477595,-2.16618e-05,4.05329e-09,-2.77406e-13,76333,-8.2527], Tmin=(816.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCJCO) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C([CH2])([O])C[CH][CH]C(10465)',
    structure = SMILES('[CH2]C([CH2])([O])C[CH][CH]C'),
    E0 = (639.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1050.26,1600,1797.85,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14997,'amu*angstrom^2'), symmetry=1, barrier=(3.44811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.519053,0.112556,-0.000176666,1.60626e-07,-5.6799e-11,77016.2,37.0278], Tmin=(100,'K'), Tmax=(860.335,'K')), NASAPolynomial(coeffs=[5.97983,0.0539244,-2.48966e-05,4.62921e-09,-3.12793e-13,76949.6,12.764], Tmin=(860.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])([CH2])[O](10466)',
    structure = SMILES('[CH2]C[CH]CC([CH2])([CH2])[O]'),
    E0 = (649.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,180,180,180,246.112,923.205,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151979,'amu*angstrom^2'), symmetry=1, barrier=(3.49429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.688313,0.113158,-0.000167857,1.44426e-07,-4.94804e-11,78324.2,36.6874], Tmin=(100,'K'), Tmax=(839.095,'K')), NASAPolynomial(coeffs=[8.95699,0.0496134,-2.28604e-05,4.2776e-09,-2.91438e-13,77324,-4.46874], Tmin=(839.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(RCCJCC) + radical(RCCJ) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C([CH2])(O)C=CCC(10467)',
    structure = SMILES('[CH2]C([CH2])(O)C=CCC'),
    E0 = (130.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.358293,0.0915424,-8.18791e-05,3.95858e-08,-7.80158e-12,15804,35.2926], Tmin=(100,'K'), Tmax=(1211.53,'K')), NASAPolynomial(coeffs=[16.0612,0.0373312,-1.47594e-05,2.65172e-09,-1.80132e-13,11825.5,-47.0724], Tmin=(1211.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(C)([O])C=CCC(10468)',
    structure = SMILES('[CH2]C(C)([O])C=CCC'),
    E0 = (145.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0213201,0.080637,-5.72114e-05,2.11663e-08,-3.22693e-12,17675.4,34.9678], Tmin=(100,'K'), Tmax=(1512.14,'K')), NASAPolynomial(coeffs=[15.9159,0.0385909,-1.5502e-05,2.77724e-09,-1.86627e-13,12868.6,-48.287], Tmin=(1512.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([CH2])(O)CC=CC(10469)',
    structure = SMILES('[CH2]C([CH2])(O)CC=CC'),
    E0 = (133.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.616283,0.105486,-0.000119584,7.58807e-08,-1.97073e-11,16215.6,32.0663], Tmin=(100,'K'), Tmax=(929.94,'K')), NASAPolynomial(coeffs=[13.9407,0.0428697,-1.85798e-05,3.46964e-09,-2.40112e-13,13508.3,-37.1051], Tmin=(929.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(C)([O])CC=CC(10470)',
    structure = SMILES('[CH2]C(C)([O])CC=CC'),
    E0 = (150.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0582455,0.0907158,-8.0466e-05,3.96414e-08,-8.20211e-12,18230.2,30.9699], Tmin=(100,'K'), Tmax=(1132.97,'K')), NASAPolynomial(coeffs=[13.0423,0.0444633,-1.92292e-05,3.60793e-09,-2.50926e-13,15261.8,-33.8681], Tmin=(1132.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = '[CH2]C([CH2])([O])C[CH]C(9431)',
    structure = SMILES('[CH2]C([CH2])([O])C[CH]C'),
    E0 = (468.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,238.649,1525.72,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.14665,'amu*angstrom^2'), symmetry=1, barrier=(3.37178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14665,'amu*angstrom^2'), symmetry=1, barrier=(3.37178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14665,'amu*angstrom^2'), symmetry=1, barrier=(3.37178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14665,'amu*angstrom^2'), symmetry=1, barrier=(3.37178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14665,'amu*angstrom^2'), symmetry=1, barrier=(3.37178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3821.38,'J/mol'), sigma=(6.83974,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.89 K, Pc=27.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0642465,0.0965236,-0.000135943,1.12228e-07,-3.73696e-11,56477.8,30.5888], Tmin=(100,'K'), Tmax=(832.788,'K')), NASAPolynomial(coeffs=[9.16048,0.0418463,-1.87814e-05,3.48587e-09,-2.36997e-13,55300.9,-10.0683], Tmin=(832.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2][C]([O])CC[CH]CC(10012)',
    structure = SMILES('[CH2][C]([O])CC[CH]CC'),
    E0 = (424.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206873,0.100346,-0.000129912,1.07885e-07,-3.73926e-11,51250.1,36.2272], Tmin=(100,'K'), Tmax=(803.41,'K')), NASAPolynomial(coeffs=[6.99963,0.0542102,-2.46281e-05,4.63171e-09,-3.18663e-13,50423.2,5.09734], Tmin=(803.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(CC)C([CH2])([CH2])[O](9471)',
    structure = SMILES('[CH2]C(CC)C([CH2])([CH2])[O]'),
    E0 = (452.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,337.624,1400.45,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3981.87,'J/mol'), sigma=(7.18225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=621.96 K, Pc=24.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.814326,0.109355,-0.000128793,8.53009e-08,-2.29654e-11,54645.4,34.3015], Tmin=(100,'K'), Tmax=(902.081,'K')), NASAPolynomial(coeffs=[14.245,0.0425766,-1.77473e-05,3.23066e-09,-2.19674e-13,51928.6,-36.7984], Tmin=(902.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)CC([CH2])([CH2])[O](10471)',
    structure = SMILES('[CH2]C(C)CC([CH2])([CH2])[O]'),
    E0 = (449.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,337.624,1400.45,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143653,'amu*angstrom^2'), symmetry=1, barrier=(3.30288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.814326,0.109355,-0.000128793,8.53009e-08,-2.29654e-11,54242.8,34.3015], Tmin=(100,'K'), Tmax=(902.081,'K')), NASAPolynomial(coeffs=[14.245,0.0425766,-1.77473e-05,3.23066e-09,-2.19674e-13,51526,-36.7984], Tmin=(902.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1(C[CH]CC)CO1(9202)',
    structure = SMILES('[CH2]C1(C[CH]CC)CO1'),
    E0 = (188.433,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3721.53,'J/mol'), sigma=(6.83365,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.29 K, Pc=26.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0413301,0.0896072,-7.63779e-05,2.3936e-08,5.19722e-12,22803.9,30.07], Tmin=(100,'K'), Tmax=(710.284,'K')), NASAPolynomial(coeffs=[11.8961,0.0414486,-1.396e-05,2.22358e-09,-1.38706e-13,20650.6,-26.3705], Tmin=(710.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJCC) + radical(CJC(C)OC)"""),
)

species(
    label = 'CC[CH]CC1([O])CC1(10472)',
    structure = SMILES('CC[CH]CC1([O])CC1'),
    E0 = (179.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.134706,0.0743714,-4.31687e-05,1.00896e-08,-2.14153e-13,21781.5,32.5315], Tmin=(100,'K'), Tmax=(1268.66,'K')), NASAPolynomial(coeffs=[14.7885,0.0396719,-1.57425e-05,2.82447e-09,-1.90883e-13,17137.6,-45.3002], Tmin=(1268.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C1([CH2])CC(CC)O1(10018)',
    structure = SMILES('[CH2]C1([CH2])CC(CC)O1'),
    E0 = (186.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12805,0.0988995,-9.5253e-05,4.98692e-08,-1.01725e-11,22604.4,30.7372], Tmin=(100,'K'), Tmax=(1331.69,'K')), NASAPolynomial(coeffs=[18.9159,0.0298359,-7.48343e-06,9.35453e-10,-4.84143e-14,18051.4,-68.7558], Tmin=(1331.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1([O])CC(CC)C1(10014)',
    structure = SMILES('[CH2]C1([O])CC(CC)C1'),
    E0 = (183.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0564212,0.070976,-1.38847e-05,-3.14154e-08,1.66098e-11,22282.1,29.1134], Tmin=(100,'K'), Tmax=(1011.48,'K')), NASAPolynomial(coeffs=[17.5958,0.0363345,-1.40008e-05,2.59702e-09,-1.84441e-13,16957.9,-64.4839], Tmin=(1011.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH]CC(6240)',
    structure = SMILES('[CH2][C]([CH2])C[CH]CC'),
    E0 = (568.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300.705,2672.4,3846.6],'cm^-1')),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735937,0.081846,-0.000108643,9.96047e-08,-3.62672e-11,68507.5,34.3431], Tmin=(100,'K'), Tmax=(876.764,'K')), NASAPolynomial(coeffs=[1.17354,0.056099,-2.39605e-05,4.31817e-09,-2.87443e-13,69343.6,37.4954], Tmin=(876.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJCC) + radical(Tertalkyl)"""),
)

species(
    label = '[CH2][C]([O])C[CH]CC(964)',
    structure = SMILES('[CH2][C]([O])C[CH]CC'),
    E0 = (448.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.454488,0.0853836,-0.000115366,9.86645e-08,-3.4579e-11,54086.8,31.6076], Tmin=(100,'K'), Tmax=(822.69,'K')), NASAPolynomial(coeffs=[6.06716,0.0458113,-2.08189e-05,3.90013e-09,-2.67059e-13,53579,8.15144], Tmin=(822.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(CJCO) + radical(RCCJCC)"""),
)

species(
    label = '[CH]CC(45)',
    structure = SMILES('[CH]CC'),
    E0 = (327.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,319.066,319.142,1779.44],'cm^-1')),
        HinderedRotor(inertia=(0.00165544,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165526,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00517,0.0155111,1.85851e-05,-3.16862e-08,1.23732e-11,39453.6,11.9344], Tmin=(100,'K'), Tmax=(982.292,'K')), NASAPolynomial(coeffs=[6.73204,0.0159276,-5.86166e-06,1.06538e-09,-7.51285e-14,37969.2,-9.8082], Tmin=(982.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])([CH2])[O](10143)',
    structure = SMILES('[CH2]C([CH2])([CH2])[O]'),
    E0 = (530.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,318.38,318.393,318.395,318.403,318.42,318.421],'cm^-1')),
        HinderedRotor(inertia=(0.0859694,'amu*angstrom^2'), symmetry=1, barrier=(6.18491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859724,'amu*angstrom^2'), symmetry=1, barrier=(6.18493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368359,'amu*angstrom^2'), symmetry=1, barrier=(26.5,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768142,0.0767143,-0.000123052,1.02056e-07,-3.2755e-11,63880.1,19.1162], Tmin=(100,'K'), Tmax=(870.983,'K')), NASAPolynomial(coeffs=[10.7052,0.020921,-9.47224e-06,1.73062e-09,-1.15193e-13,62534.4,-25.24], Tmin=(870.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C[CH2](6)',
    structure = SMILES('C[CH2]'),
    E0 = (108.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,474.132,1048.55,2319.88,2320.55,2321.73],'cm^-1')),
        HinderedRotor(inertia=(0.000749852,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82183,-0.00343357,5.09256e-05,-6.2021e-08,2.37073e-11,13066,7.61644], Tmin=(100,'K'), Tmax=(900.314,'K')), NASAPolynomial(coeffs=[5.15622,0.00943121,-1.81945e-06,2.21194e-10,-1.4348e-14,12064.1,-2.91103], Tmin=(900.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ)"""),
)

species(
    label = '[CH]CC([CH2])([CH2])[O](10283)',
    structure = SMILES('[CH]CC([CH2])([CH2])[O]'),
    E0 = (745.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,891.781,1333.73,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159229,'amu*angstrom^2'), symmetry=1, barrier=(3.66098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159229,'amu*angstrom^2'), symmetry=1, barrier=(3.66098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159229,'amu*angstrom^2'), symmetry=1, barrier=(3.66098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159229,'amu*angstrom^2'), symmetry=1, barrier=(3.66098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471788,0.0826618,-0.000113906,8.17394e-08,-2.24654e-11,89842.2,24.7037], Tmin=(100,'K'), Tmax=(717.367,'K')), NASAPolynomial(coeffs=[12.1227,0.0256988,-1.15291e-05,2.1477e-09,-1.469e-13,87964.7,-29.0703], Tmin=(717.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(745.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CJC(C)2O) + radical(CJC(C)2O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])([O])C[C]CC(10473)',
    structure = SMILES('[CH2]C([CH2])([O])C[C]CC'),
    E0 = (698.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,859.413,1360.84,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157139,'amu*angstrom^2'), symmetry=1, barrier=(3.61293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.849245,0.112632,-0.000143563,1.01727e-07,-2.9268e-11,84168.8,32.8347], Tmin=(100,'K'), Tmax=(845.243,'K')), NASAPolynomial(coeffs=[13.9646,0.0425303,-1.91642e-05,3.61434e-09,-2.50365e-13,81664.4,-36.1437], Tmin=(845.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(698.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CJC(C)2O) + radical(CCJ2_triplet) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH]C([CH2])([O])C[CH]CC(10474)',
    structure = SMILES('[CH]C([CH2])([O])C[CH]CC'),
    E0 = (680.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.426132,0.104367,-0.000138127,1.1042e-07,-3.64994e-11,82040.5,35.8218], Tmin=(100,'K'), Tmax=(794.802,'K')), NASAPolynomial(coeffs=[9.58995,0.0483777,-2.19253e-05,4.1166e-09,-2.83011e-13,80624.7,-9.0906], Tmin=(794.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)2O) + radical(CCJ2_triplet) + radical(CC(C)2OJ) + radical(RCCJCC)"""),
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
    E0 = (444.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (576.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (579.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (444.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (485.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (444.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (535.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (556.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (603.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (603.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (519.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (562.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (596.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (495.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (491.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (538.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (497.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (792.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (856.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (809.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (850.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (861.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (508.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (508.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (469.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (469.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (887.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (601.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (612.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (644.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (450.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (452.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (452.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (452.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (975.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (864.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (892.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (888.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (910.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (892.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['C=CCC(36)', '[CH2]C(=C)[O](4273)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2]C([CH2])([O])C=CCC(10454)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C([CH2])([O])CC=CC(10455)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C[CH]CC(950)', 'CH2(T)(28)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(9.66262,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [CO_O;CH2_triplet]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from -5.8 to 9.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(T)(63)', '[CH2]C(=C)C[CH]CC(6234)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]CC(39)', '[CH2]C(=C)[O](4273)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(100.98,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 100.1 to 101.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CCC(36)', '[CH2][C]([CH2])[O](10271)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0172287,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH3](11)', '[CH2]C([CH2])([O])CC=C(10187)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])([O])[CH]CCC(10423)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])([O])CC[CH]C(10456)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])(O)[CH][CH]CC(10457)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C(C)([O])[CH][CH]CC(10458)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CCCC([CH2])([CH2])[O](10459)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])(O)C[CH][CH]C(10460)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C(C)([O])C[CH][CH]C(10461)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(524000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C[CH]CC([CH2])([CH2])O(10462)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_2H] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C[CH]CC([CH2])(C)[O](10463)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(10873.3,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R6HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]CC(39)', '[CH2][C]([CH2])[O](10271)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2]C([CH2])([O])[CH][CH]CC(10464)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH3](11)', '[CH2][CH]CC([CH2])([CH2])[O](9961)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.66881e+08,'m^3/(mol*s)'), n=-0.401267, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -6.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[CH2]C([CH2])([O])C[CH][CH]C(10465)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH2]C[CH]CC([CH2])([CH2])[O](10466)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])(O)C=CCC(10467)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C(C)([O])C=CCC(10468)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C([CH2])(O)CC=CC(10469)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C(C)([O])CC=CC(10470)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.50442e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(14)', '[CH2]C([CH2])([O])C[CH]C(9431)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2][C]([O])CC[CH]CC(10012)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CC)C([CH2])([CH2])[O](9471)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)CC([CH2])([CH2])[O](10471)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C1(C[CH]CC)CO1(9202)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.18842e+14,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['CC[CH]CC1([O])CC1(10472)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C1([CH2])CC(CC)O1(10018)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    products = ['[CH2]C1([O])CC(CC)C1(10014)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', '[CH2][C]([CH2])C[CH]CC(6240)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(T)(28)', '[CH2][C]([O])C[CH]CC(964)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]CC(45)', '[CH2]C([CH2])([CH2])[O](10143)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.44562e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[CH2](6)', '[CH]CC([CH2])([CH2])[O](10283)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(8)', '[CH2]C([CH2])([O])C[C]CC(10473)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(8)', '[CH]C([CH2])([O])C[CH]CC(10474)'],
    products = ['[CH2]C([CH2])([O])C[CH]CC(9519)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2624',
    isomers = [
        '[CH2]C([CH2])([O])C[CH]CC(9519)',
    ],
    reactants = [
        ('C=CCC(36)', '[CH2]C(=C)[O](4273)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2624',
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


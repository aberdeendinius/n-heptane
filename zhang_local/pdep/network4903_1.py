species(
    label = 'C#CC([O])C([CH2])=C(22602)',
    structure = SMILES('C#CC([O])C([CH2])=C'),
    E0 = (419.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,381.987,382.066],'cm^-1')),
        HinderedRotor(inertia=(0.270748,'amu*angstrom^2'), symmetry=1, barrier=(28.0456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270732,'amu*angstrom^2'), symmetry=1, barrier=(28.0431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449297,'amu*angstrom^2'), symmetry=1, barrier=(46.5425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83795,0.0613809,-4.54481e-05,9.35844e-09,2.54232e-12,50596,25.0188], Tmin=(100,'K'), Tmax=(975.246,'K')), NASAPolynomial(coeffs=[15.5972,0.0184084,-6.3662e-06,1.10825e-09,-7.63922e-14,46882,-50.098], Tmin=(975.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=C(6077)',
    structure = SMILES('C=C=C'),
    E0 = (182.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36163,0.00777114,2.52438e-05,-3.61894e-08,1.38594e-11,22005.7,7.1245], Tmin=(100,'K'), Tmax=(966.19,'K')), NASAPolynomial(coeffs=[6.46487,0.0106812,-3.73734e-06,6.87011e-10,-4.98581e-14,20670.6,-11.5463], Tmin=(966.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C#CC1OC1([CH2])[CH2](24397)',
    structure = SMILES('C#CC1OC1([CH2])[CH2]'),
    E0 = (509.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.472439,0.0882854,-0.000113939,7.1391e-08,-1.63314e-11,61390.7,24.7541], Tmin=(100,'K'), Tmax=(1286.22,'K')), NASAPolynomial(coeffs=[18.6208,0.00882971,2.13798e-06,-9.09629e-10,8.04478e-14,58139.9,-65.7097], Tmin=(1286.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C1OC1C([CH2])=C(24401)',
    structure = SMILES('[CH]=C1OC1C([CH2])=C'),
    E0 = (441.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875227,0.0488895,1.19047e-05,-6.57184e-08,3.27887e-11,53261.5,21.7229], Tmin=(100,'K'), Tmax=(937.166,'K')), NASAPolynomial(coeffs=[22.7007,0.006016,4.72722e-08,-3.3179e-11,-5.79263e-15,46962.6,-93.9368], Tmin=(937.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C1CC(=C)C1[O](24287)',
    structure = SMILES('[CH]=C1CC(=C)C1[O]'),
    E0 = (495.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72839,0.0308545,4.54732e-05,-8.50756e-08,3.54676e-11,59725.2,23.7722], Tmin=(100,'K'), Tmax=(964.524,'K')), NASAPolynomial(coeffs=[16.7555,0.0160145,-5.28666e-06,1.04533e-09,-8.29604e-14,54617.9,-59.631], Tmin=(964.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Cds_P)"""),
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
    label = 'C#CC(=O)C([CH2])=C(24402)',
    structure = SMILES('C#CC(=O)C([CH2])=C'),
    E0 = (255.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2175,525,375,552.5,462.5,1710,311.317],'cm^-1')),
        HinderedRotor(inertia=(0.00173916,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6222,'amu*angstrom^2'), symmetry=1, barrier=(42.8177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843746,'amu*angstrom^2'), symmetry=1, barrier=(58.0193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52453,0.0560387,-5.39551e-05,2.90545e-08,-6.54185e-12,30815.6,20.0933], Tmin=(100,'K'), Tmax=(1051.34,'K')), NASAPolynomial(coeffs=[9.33184,0.0263353,-1.15771e-05,2.18291e-09,-1.52195e-13,29174,-17.9636], Tmin=(1051.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2][C]=C(6078)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.243012,'amu*angstrom^2'), symmetry=1, barrier=(30.4931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28772,0.0102496,1.79964e-05,-2.86352e-08,1.11956e-11,47593.9,10.4766], Tmin=(100,'K'), Tmax=(969.996,'K')), NASAPolynomial(coeffs=[6.37267,0.0109726,-3.91218e-06,7.11399e-10,-5.07602e-14,46362.9,-7.57277], Tmin=(969.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH2]C(=C)C=O(17490)',
    structure = SMILES('[CH2]C(=C)C=O'),
    E0 = (31.6694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0868503,'amu*angstrom^2'), symmetry=1, barrier=(19.2321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836313,'amu*angstrom^2'), symmetry=1, barrier=(19.2285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12629,0.0353729,-1.91494e-05,2.12252e-09,8.03052e-13,3881.16,15.3367], Tmin=(100,'K'), Tmax=(1283.06,'K')), NASAPolynomial(coeffs=[10.8189,0.0169683,-7.79819e-06,1.50631e-09,-1.06015e-13,934.834,-31.5557], Tmin=(1283.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.6694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C=O)CJ)"""),
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
    label = 'C#CC(O)=C([CH2])[CH2](24403)',
    structure = SMILES('C#CC(O)=C([CH2])[CH2]'),
    E0 = (224.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564542,0.0711446,-7.14875e-05,3.18341e-08,-3.49117e-12,27137,20.3288], Tmin=(100,'K'), Tmax=(888.237,'K')), NASAPolynomial(coeffs=[15.922,0.0168606,-4.93656e-06,7.38423e-10,-4.57041e-14,23822,-55.2451], Tmin=(888.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = 'C#CC([O])=C([CH2])C(24404)',
    structure = SMILES('C#CC([O])=C([CH2])C'),
    E0 = (240.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.836272,0.0685508,-7.53738e-05,4.48658e-08,-1.06783e-11,29095,21.2032], Tmin=(100,'K'), Tmax=(1024.29,'K')), NASAPolynomial(coeffs=[12.4781,0.0230877,-8.79631e-06,1.53328e-09,-1.02043e-13,26710.1,-35.2413], Tmin=(1024.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=C(C)C([O])C#C(24405)',
    structure = SMILES('[CH]=C(C)C([O])C#C'),
    E0 = (515.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1380,1390,370,380,2900,435,2175,525,193.164],'cm^-1')),
        HinderedRotor(inertia=(0.660871,'amu*angstrom^2'), symmetry=1, barrier=(17.4985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40099,'amu*angstrom^2'), symmetry=1, barrier=(37.0985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66082,'amu*angstrom^2'), symmetry=1, barrier=(17.4983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827513,0.0666449,-6.79795e-05,3.65501e-08,-7.8295e-12,62089.4,25.8728], Tmin=(100,'K'), Tmax=(1134.71,'K')), NASAPolynomial(coeffs=[13.6438,0.0214655,-8.25555e-06,1.46083e-09,-9.85585e-14,59180.9,-37.578], Tmin=(1134.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C(O)C#C(24406)',
    structure = SMILES('[CH]C(=C)C(O)C#C'),
    E0 = (408.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611596,0.0685209,-6.19107e-05,2.9979e-08,-5.85294e-12,49257.8,27.5455], Tmin=(100,'K'), Tmax=(1234.23,'K')), NASAPolynomial(coeffs=[13.9491,0.0252962,-9.37901e-06,1.60445e-09,-1.05599e-13,45965.4,-39.6073], Tmin=(1234.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#CC(O)C([CH2])=C(24407)',
    structure = SMILES('[C]#CC(O)C([CH2])=C'),
    E0 = (526.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2175,525,1380,1390,370,380,2900,435,180,487.151],'cm^-1')),
        HinderedRotor(inertia=(0.0702872,'amu*angstrom^2'), symmetry=1, barrier=(11.8368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425529,'amu*angstrom^2'), symmetry=1, barrier=(71.6622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702875,'amu*angstrom^2'), symmetry=1, barrier=(11.8368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425527,'amu*angstrom^2'), symmetry=1, barrier=(71.6622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735358,0.0684747,-7.23743e-05,3.92692e-08,-8.00613e-12,63438.4,26.3214], Tmin=(100,'K'), Tmax=(942.254,'K')), NASAPolynomial(coeffs=[13.8965,0.0199335,-6.76856e-06,1.10726e-09,-7.13371e-14,60632.8,-38.1172], Tmin=(942.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC([O])C(=C)C(24408)',
    structure = SMILES('[C]#CC([O])C(=C)C'),
    E0 = (605.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,350,440,435,1725,264.14,264.142,264.144],'cm^-1')),
        HinderedRotor(inertia=(0.00241616,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322578,'amu*angstrom^2'), symmetry=1, barrier=(15.9715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40753,'amu*angstrom^2'), symmetry=1, barrier=(69.6856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.977515,0.0660179,-7.1581e-05,4.30417e-08,-1.04431e-11,72911.9,25.6752], Tmin=(100,'K'), Tmax=(1001.76,'K')), NASAPolynomial(coeffs=[11.3847,0.0244615,-9.35512e-06,1.63006e-09,-1.08283e-13,70826.9,-24.5518], Tmin=(1001.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([CH2])=C[O](8501)',
    structure = SMILES('[CH2]C([CH2])=C[O]'),
    E0 = (202.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(1.66675,'amu*angstrom^2'), symmetry=1, barrier=(38.3218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66433,'amu*angstrom^2'), symmetry=1, barrier=(38.2662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86414,0.0317844,2.03304e-05,-6.01192e-08,2.85015e-11,24474.6,17.4197], Tmin=(100,'K'), Tmax=(926.323,'K')), NASAPolynomial(coeffs=[17.1232,0.00501783,3.19565e-07,-1.22369e-10,3.89574e-15,19969.1,-64.0889], Tmin=(926.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([O])=C([CH2])[CH2](24409)',
    structure = SMILES('C#CC([O])=C([CH2])[CH2]'),
    E0 = (362.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,325,375,415,465,420,450,1700,1750,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,286.903],'cm^-1')),
        HinderedRotor(inertia=(0.00214449,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00315,'amu*angstrom^2'), symmetry=1, barrier=(57.898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01491,'amu*angstrom^2'), symmetry=1, barrier=(57.946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00083,0.0662531,-7.53238e-05,4.33695e-08,-8.7386e-12,43691.1,20.2698], Tmin=(100,'K'), Tmax=(804.28,'K')), NASAPolynomial(coeffs=[12.0847,0.0204018,-7.10426e-06,1.15764e-09,-7.34844e-14,41608.3,-32.6537], Tmin=(804.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]C(=C)C([O])C#C(22616)',
    structure = SMILES('[CH]C(=C)C([O])C#C'),
    E0 = (638.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,350,440,435,1725,290.584,290.77,290.899,290.964,291.327],'cm^-1')),
        HinderedRotor(inertia=(0.838736,'amu*angstrom^2'), symmetry=1, barrier=(50.1572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.839241,'amu*angstrom^2'), symmetry=1, barrier=(50.1609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835541,'amu*angstrom^2'), symmetry=1, barrier=(50.1597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725246,0.0650824,-5.7831e-05,2.71461e-08,-5.11573e-12,76960.4,26.9472], Tmin=(100,'K'), Tmax=(1279.3,'K')), NASAPolynomial(coeffs=[14.1709,0.0230419,-8.53778e-06,1.45846e-09,-9.58735e-14,73520.3,-41.232], Tmin=(1279.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])=C(24410)',
    structure = SMILES('[C]#CC([O])C([CH2])=C'),
    E0 = (756.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,350,440,435,1725,392.591,400.99,847.619],'cm^-1')),
        HinderedRotor(inertia=(0.644594,'amu*angstrom^2'), symmetry=1, barrier=(71.2588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151124,'amu*angstrom^2'), symmetry=1, barrier=(17.2939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.729361,'amu*angstrom^2'), symmetry=1, barrier=(71.2612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884332,0.0646496,-6.70837e-05,3.5069e-08,-6.77308e-12,91139.5,25.5945], Tmin=(100,'K'), Tmax=(957.72,'K')), NASAPolynomial(coeffs=[13.8742,0.0180673,-6.14025e-06,1.00971e-09,-6.55156e-14,88299.6,-38.3498], Tmin=(957.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(756.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Acetyl)"""),
)

species(
    label = 'C#CC([O])[C]1CC1(24411)',
    structure = SMILES('C#CC([O])[C]1CC1'),
    E0 = (465.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61875,0.0340419,3.88937e-05,-8.25048e-08,3.63412e-11,56095.4,23.2744], Tmin=(100,'K'), Tmax=(935.773,'K')), NASAPolynomial(coeffs=[17.3109,0.0135448,-2.9157e-06,4.74497e-10,-3.84629e-14,51119.1,-62.2867], Tmin=(935.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]=[C]C1OCC1=C(24323)',
    structure = SMILES('[CH]=[C]C1OCC1=C'),
    E0 = (552.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69665,0.0385078,8.88961e-06,-3.67534e-08,1.57456e-11,66535.1,24.4095], Tmin=(100,'K'), Tmax=(1026.48,'K')), NASAPolynomial(coeffs=[12.8067,0.0232167,-9.68532e-06,1.88639e-09,-1.37749e-13,62779,-36.6662], Tmin=(1026.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C)C1[C]=CO1(24378)',
    structure = SMILES('[CH2]C(=C)C1[C]=CO1'),
    E0 = (413.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05157,0.0408095,4.10421e-05,-9.86001e-08,4.49665e-11,49804.8,21.7334], Tmin=(100,'K'), Tmax=(935.772,'K')), NASAPolynomial(coeffs=[23.858,0.00453805,1.05787e-06,-2.07273e-10,3.62719e-15,42856.3,-101.101], Tmin=(935.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1CC=[C]C1[O](24311)',
    structure = SMILES('C=C1CC=[C]C1[O]'),
    E0 = (425.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03134,0.0300903,2.89145e-05,-5.44892e-08,2.13033e-11,51246.8,23.0441], Tmin=(100,'K'), Tmax=(1017.62,'K')), NASAPolynomial(coeffs=[11.5003,0.0247798,-1.02935e-05,2.01123e-09,-1.47485e-13,47667.5,-30.9213], Tmin=(1017.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC(=O)C(=C)C(24412)',
    structure = SMILES('C#CC(=O)C(=C)C'),
    E0 = (99.1339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80879,0.0541731,-4.42196e-05,9.72334e-09,8.53464e-12,11996.3,20.3174], Tmin=(100,'K'), Tmax=(564.128,'K')), NASAPolynomial(coeffs=[5.64918,0.0353221,-1.63766e-05,3.15061e-09,-2.21545e-13,11429.7,2.80644], Tmin=(564.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.1339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])[CH2](23626)',
    structure = SMILES('[CH]=C=C([O])C([CH2])[CH2]'),
    E0 = (580.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,350,440,435,1725,259.688],'cm^-1')),
        HinderedRotor(inertia=(1.59674,'amu*angstrom^2'), symmetry=1, barrier=(76.6401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159914,'amu*angstrom^2'), symmetry=1, barrier=(7.52698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0117413,'amu*angstrom^2'), symmetry=1, barrier=(76.6779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41285,0.0716126,-8.23635e-05,4.91474e-08,-1.10865e-11,69916.5,30.8561], Tmin=(100,'K'), Tmax=(1250.7,'K')), NASAPolynomial(coeffs=[14.9478,0.0150468,-2.43329e-06,9.77769e-11,6.00553e-15,67069.1,-39.3663], Tmin=(1250.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Isobutyl) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])[CH2](23627)',
    structure = SMILES('[C]#CC([O])C([CH2])[CH2]'),
    E0 = (897.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,394.032,394.032,2232.91],'cm^-1')),
        HinderedRotor(inertia=(0.645897,'amu*angstrom^2'), symmetry=1, barrier=(71.163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645896,'amu*angstrom^2'), symmetry=1, barrier=(71.163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00264042,'amu*angstrom^2'), symmetry=1, barrier=(9.34204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645895,'amu*angstrom^2'), symmetry=1, barrier=(71.163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.555856,0.0715521,-8.45548e-05,5.31556e-08,-1.27831e-11,108095,30.9446], Tmin=(100,'K'), Tmax=(1160.04,'K')), NASAPolynomial(coeffs=[13.5206,0.0180927,-4.10798e-06,4.17457e-10,-1.54152e-14,105676,-30.9882], Tmin=(1160.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(897.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(CC(C)OJ) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = '[C]#CC([O])[C]([CH2])C(24413)',
    structure = SMILES('[C]#CC([O])[C]([CH2])C'),
    E0 = (845.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2175,525,360,370,350,456.207,456.256,1993.25],'cm^-1')),
        HinderedRotor(inertia=(0.0775047,'amu*angstrom^2'), symmetry=1, barrier=(11.4464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025475,'amu*angstrom^2'), symmetry=1, barrier=(71.8238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486228,'amu*angstrom^2'), symmetry=1, barrier=(71.827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486401,'amu*angstrom^2'), symmetry=1, barrier=(71.8213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736337,0.06401,-6.31442e-05,3.36211e-08,-7.00656e-12,101776,29.4081], Tmin=(100,'K'), Tmax=(1274.91,'K')), NASAPolynomial(coeffs=[13.726,0.0190056,-5.19445e-06,7.04065e-10,-3.91185e-14,98809.1,-35.0601], Tmin=(1274.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(Acetyl) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C#CC([O])C[C]=C(22605)',
    structure = SMILES('C#CC([O])C[C]=C'),
    E0 = (525.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,2175,525,390.21,390.238,390.243],'cm^-1')),
        HinderedRotor(inertia=(0.335851,'amu*angstrom^2'), symmetry=1, barrier=(36.2732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150774,'amu*angstrom^2'), symmetry=1, barrier=(16.2912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150759,'amu*angstrom^2'), symmetry=1, barrier=(16.2932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3808.29,'J/mol'), sigma=(6.33679,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.85 K, Pc=33.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77991,0.0635805,-6.01783e-05,2.96196e-08,-5.75889e-12,63270.5,27.7658], Tmin=(100,'K'), Tmax=(1253.96,'K')), NASAPolynomial(coeffs=[14.8428,0.0187217,-6.51823e-06,1.09153e-09,-7.13415e-14,59743.6,-43.2622], Tmin=(1253.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1OCC1=C(24336)',
    structure = SMILES('C#CC1OCC1=C'),
    E0 = (233.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79344,0.0331401,3.15149e-05,-6.6268e-08,2.81283e-11,28169.6,20.7991], Tmin=(100,'K'), Tmax=(961.933,'K')), NASAPolynomial(coeffs=[14.2122,0.0195261,-6.55387e-06,1.21186e-09,-8.99837e-14,24021,-47.7774], Tmin=(961.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2][C]1CC=[C]C1[O](24414)',
    structure = SMILES('[CH2][C]1CC=[C]C1[O]'),
    E0 = (661.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01318,0.0260913,5.18125e-05,-8.78893e-08,3.60399e-11,79694.8,24.5932], Tmin=(100,'K'), Tmax=(952.634,'K')), NASAPolynomial(coeffs=[14.372,0.0183971,-5.66763e-06,1.03978e-09,-7.90591e-14,75334.5,-44.9579], Tmin=(952.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CCJ(C)CO) + radical(CC(C)OJ) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = 'C#CC=C([CH2])[CH2](19657)',
    structure = SMILES('C#CC=C([CH2])[CH2]'),
    E0 = (438.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.142552,'amu*angstrom^2'), symmetry=1, barrier=(71.2344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.10175,'amu*angstrom^2'), symmetry=1, barrier=(71.3153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.10514,'amu*angstrom^2'), symmetry=1, barrier=(71.3934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39256,0.0509418,-4.21539e-05,1.91619e-08,-3.53282e-12,52833.8,18.1133], Tmin=(100,'K'), Tmax=(1303.18,'K')), NASAPolynomial(coeffs=[11.377,0.0202952,-6.87877e-06,1.11622e-09,-7.09586e-14,50231.5,-32.7], Tmin=(1303.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
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
    label = 'C#CC([O])[C]=C(23414)',
    structure = SMILES('C#CC([O])[C]=C'),
    E0 = (545.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,1685,370,255.784,258.255],'cm^-1')),
        HinderedRotor(inertia=(0.729244,'amu*angstrom^2'), symmetry=1, barrier=(33.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519157,'amu*angstrom^2'), symmetry=1, barrier=(24.3049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63247,0.0486645,-4.88996e-05,2.57615e-08,-5.38025e-12,65644.8,22.2968], Tmin=(100,'K'), Tmax=(1166.95,'K')), NASAPolynomial(coeffs=[11.494,0.0148618,-5.44949e-06,9.38862e-10,-6.23884e-14,63343.2,-26.8018], Tmin=(1166.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(CC(C)OJ)"""),
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
    E0 = (419.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (512.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (504.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (495.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (488.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (499.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (641.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (484.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (545.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (622.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (644.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (839.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (610.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (645.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (665.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (754.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (574.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (851.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (968.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (650.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (552.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (661.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (478.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (483.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (603.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (922.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (870.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (695.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (427.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (812.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (845.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (960.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C=C=C(6077)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC1OC1([CH2])[CH2](24397)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.33596e+11,'s^-1'), n=-0.0500183, Ea=(92.3582,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['[CH]=C1OC1C([CH2])=C(24401)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['[CH]=C1CC(=C)C1[O](24287)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(354414,'s^-1'), n=1.88643, Ea=(76.0883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 72.1 to 76.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#CC(=O)C([CH2])=C(24402)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][C]=C(6078)', 'C#CC=O(21959)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.99e+06,'cm^3/(mol*s)'), n=2.12, Ea=(20.0832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CtH_O;YJ] for rate rule [CO-CtH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]#C(5143)', '[CH2]C(=C)C=O(17490)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C=C(6077)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC(O)=C([CH2])[CH2](24403)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC([O])=C([CH2])C(24404)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C)C([O])C#C(24405)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['[CH]C(=C)C(O)C#C(24406)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.84141e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC(O)C([CH2])=C(24407)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CC([O])C(=C)C(24408)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]=C(6078)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#C(5143)', '[CH2]C([CH2])=C[O](8501)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', 'C#CC([O])=C([CH2])[CH2](24409)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]C(=C)C([O])C#C(22616)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[C]#CC([O])C([CH2])=C(24410)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC([O])[C]1CC1(24411)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['[CH]=[C]C1OCC1=C(24323)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.06838e+08,'s^-1'), n=1.06803, Ea=(132.755,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 129.5 to 132.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['[CH2]C(=C)C1[C]=CO1(24378)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C=C1CC=[C]C1[O](24311)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.94e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC(=O)C(=C)C(24412)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C=C([O])C([CH2])[CH2](23626)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[C]#CC([O])C([CH2])[CH2](23627)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[C]#CC([O])[C]([CH2])C(24413)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([O])C[C]=C(22605)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([O])C([CH2])=C(22602)'],
    products = ['C#CC1OCC1=C(24336)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]1CC=[C]C1[O](24414)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(T)(63)', 'C#CC=C([CH2])[CH2](19657)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(28)', 'C#CC([O])[C]=C(23414)'],
    products = ['C#CC([O])C([CH2])=C(22602)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

network(
    label = '4903',
    isomers = [
        'C#CC([O])C([CH2])=C(22602)',
    ],
    reactants = [
        ('C=C=C(6077)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4903',
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


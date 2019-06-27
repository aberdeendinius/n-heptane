species(
    label = '[CH2][C]=CC[C]([O])[O](18581)',
    structure = SMILES('[CH2][C]=CC[C]([O])[O]'),
    E0 = (633.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,360,370,350,180,1240.91,1241.17,1241.21],'cm^-1')),
        HinderedRotor(inertia=(0.168189,'amu*angstrom^2'), symmetry=1, barrier=(3.86699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168191,'amu*angstrom^2'), symmetry=1, barrier=(3.86705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168292,'amu*angstrom^2'), symmetry=1, barrier=(3.86937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42851,0.0668148,-0.000110334,1.10566e-07,-4.2961e-11,76262,29.4431], Tmin=(100,'K'), Tmax=(822.245,'K')), NASAPolynomial(coeffs=[2.12329,0.0405539,-2.06857e-05,4.03706e-09,-2.81357e-13,76921.3,30.9306], Tmin=(822.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S) + radical(Allyl_P) + radical(Cs_P)"""),
)

species(
    label = 'C=C([O])[O](1172)',
    structure = SMILES('C=C([O])[O]'),
    E0 = (-40.8548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,219.703,219.72],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3790.78,'J/mol'), sigma=(6.02099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.11 K, Pc=39.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85742,0.0201906,-6.55878e-06,-8.62627e-09,5.39429e-12,-4868.13,11.9276], Tmin=(100,'K'), Tmax=(974.512,'K')), NASAPolynomial(coeffs=[9.20068,0.00568109,-1.96826e-06,3.71357e-10,-2.78204e-14,-6651.8,-21.3196], Tmin=(974.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.8548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[CH2][C]=CC=C([O])[O](20892)',
    structure = SMILES('[CH2][C]=CC=C([O])[O]'),
    E0 = (330.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,319.744,319.845,319.847],'cm^-1')),
        HinderedRotor(inertia=(0.49416,'amu*angstrom^2'), symmetry=1, barrier=(35.8355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493366,'amu*angstrom^2'), symmetry=1, barrier=(35.8297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12836,0.0563532,-4.79461e-05,1.58779e-08,-3.29511e-13,39894.2,24.6061], Tmin=(100,'K'), Tmax=(987.529,'K')), NASAPolynomial(coeffs=[14.7884,0.0145983,-5.14281e-06,9.02251e-10,-6.23552e-14,36534.3,-44.4757], Tmin=(987.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C#CC[C]([O])[O](20893)',
    structure = SMILES('[CH2]C#CC[C]([O])[O]'),
    E0 = (560.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2100,2250,500,550,360,370,350,180,1043.1,1043.17,1043.25],'cm^-1')),
        HinderedRotor(inertia=(0.151404,'amu*angstrom^2'), symmetry=1, barrier=(3.48107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15143,'amu*angstrom^2'), symmetry=1, barrier=(3.48167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00450772,'amu*angstrom^2'), symmetry=1, barrier=(3.48167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72565,0.0601802,-0.000101328,1.01837e-07,-3.89088e-11,67473.3,28.1384], Tmin=(100,'K'), Tmax=(850.126,'K')), NASAPolynomial(coeffs=[1.57779,0.0369416,-1.8094e-05,3.44776e-09,-2.36137e-13,68363.3,33.9145], Tmin=(850.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCOJ) + radical(Cs_P) + radical(Propargyl) + radical(CCOJ)"""),
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
    label = '[CH2][C]([O])[O](2304)',
    structure = SMILES('[CH2][C]([O])[O]'),
    E0 = (411.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,1740.05,1740.26],'cm^-1')),
        HinderedRotor(inertia=(0.00349141,'amu*angstrom^2'), symmetry=1, barrier=(7.50228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01691,0.0315826,-7.13533e-05,8.28598e-08,-3.37644e-11,49511,16.8536], Tmin=(100,'K'), Tmax=(857.006,'K')), NASAPolynomial(coeffs=[-0.765291,0.0236786,-1.27872e-05,2.50409e-09,-1.7282e-13,51097.8,39.9925], Tmin=(857.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH2][C]=C[CH]C([O])[O](20894)',
    structure = SMILES('[CH2][C]=C[CH]C([O])[O]'),
    E0 = (545.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20558,0.0696491,-0.000106309,1.00607e-07,-3.84026e-11,65648.5,26.755], Tmin=(100,'K'), Tmax=(791.029,'K')), NASAPolynomial(coeffs=[4.26705,0.0391677,-2.00638e-05,3.94823e-09,-2.77665e-13,65633.4,15.6693], Tmin=(791.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S) + radical(Allyl_P) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C=[C]C[C]([O])[O](20895)',
    structure = SMILES('[CH2]C=[C]C[C]([O])[O]'),
    E0 = (633.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,360,370,350,180,1239.52,1240.41,1243.53],'cm^-1')),
        HinderedRotor(inertia=(0.166341,'amu*angstrom^2'), symmetry=1, barrier=(3.8245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169269,'amu*angstrom^2'), symmetry=1, barrier=(3.89183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16914,'amu*angstrom^2'), symmetry=1, barrier=(3.88885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42854,0.0668143,-0.000110332,1.10563e-07,-4.29596e-11,76262,29.443], Tmin=(100,'K'), Tmax=(822.251,'K')), NASAPolynomial(coeffs=[2.12323,0.040554,-2.06858e-05,4.03708e-09,-2.81359e-13,76921.3,30.9309], Tmin=(822.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=[C]CC([O])[O](20896)',
    structure = SMILES('[CH2][C]=[C]CC([O])[O]'),
    E0 = (665.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1229.09,1229.86,1229.89],'cm^-1')),
        HinderedRotor(inertia=(0.22569,'amu*angstrom^2'), symmetry=1, barrier=(5.18905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225499,'amu*angstrom^2'), symmetry=1, barrier=(5.18467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225747,'amu*angstrom^2'), symmetry=1, barrier=(5.19037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34828,0.0707982,-0.000125062,1.28021e-07,-4.96199e-11,80183.1,29.5168], Tmin=(100,'K'), Tmax=(838.044,'K')), NASAPolynomial(coeffs=[1.22836,0.0421941,-2.16413e-05,4.207e-09,-2.91478e-13,81227.8,36.1869], Tmin=(838.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH][CH]C=C([O])[O](20897)',
    structure = SMILES('[CH2][CH][CH]C=C([O])[O]'),
    E0 = (417.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63924,0.0546299,-5.17011e-05,2.78172e-08,-6.36933e-12,50296.3,27.3801], Tmin=(100,'K'), Tmax=(1023.77,'K')), NASAPolynomial(coeffs=[8.44802,0.0280271,-1.27233e-05,2.43535e-09,-1.71177e-13,48902.2,-5.62831], Tmin=(1023.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(RCCJ) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJC) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]=C[CH][C]([O])O(20898)',
    structure = SMILES('[CH2][C]=C[CH][C]([O])O'),
    E0 = (524.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01093,0.0713919,-9.73291e-05,7.67202e-08,-2.50678e-11,63197.6,27.9866], Tmin=(100,'K'), Tmax=(741.373,'K')), NASAPolynomial(coeffs=[8.5177,0.0308877,-1.53734e-05,3.01882e-09,-2.13435e-13,62084.6,-5.98235], Tmin=(741.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Cs_P) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C[C]=[C]C[C]([O])[O](20899)',
    structure = SMILES('C[C]=[C]C[C]([O])[O]'),
    E0 = (719.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,2324.57,2327.98],'cm^-1')),
        HinderedRotor(inertia=(0.160602,'amu*angstrom^2'), symmetry=1, barrier=(3.69256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160079,'amu*angstrom^2'), symmetry=1, barrier=(3.68053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160971,'amu*angstrom^2'), symmetry=1, barrier=(3.70104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37543,0.0741606,-0.00014461,1.54462e-07,-6.04756e-11,86642.7,29.8261], Tmin=(100,'K'), Tmax=(853.099,'K')), NASAPolynomial(coeffs=[-1.17885,0.045793,-2.37945e-05,4.61473e-09,-3.17783e-13,88546.6,50.3476], Tmin=(853.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C[C]([O])O(20900)',
    structure = SMILES('[CH2][C]=[C]C[C]([O])O'),
    E0 = (645.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1670,1700,300,440,360,370,350,180,919.794,921.361],'cm^-1')),
        HinderedRotor(inertia=(0.162298,'amu*angstrom^2'), symmetry=1, barrier=(3.73154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163339,'amu*angstrom^2'), symmetry=1, barrier=(3.75548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162585,'amu*angstrom^2'), symmetry=1, barrier=(3.73816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163499,'amu*angstrom^2'), symmetry=1, barrier=(3.75916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04397,0.0740045,-0.000122092,1.13416e-07,-4.10205e-11,77736.9,31.1314], Tmin=(100,'K'), Tmax=(834.706,'K')), NASAPolynomial(coeffs=[5.71944,0.0334731,-1.66826e-05,3.21175e-09,-2.21629e-13,77587.8,13.2021], Tmin=(834.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cs_P) + radical(Cds_S) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C][CH]C=C([O])[O](20901)',
    structure = SMILES('C[C][CH]C=C([O])[O]'),
    E0 = (466.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.993737,0.0610157,-5.64239e-05,2.62928e-08,-4.86714e-12,56161.9,24.3812], Tmin=(100,'K'), Tmax=(1301.93,'K')), NASAPolynomial(coeffs=[14.6827,0.0189582,-7.96792e-06,1.48037e-09,-1.0259e-13,52597.5,-45.2718], Tmin=(1301.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(CCJ2_triplet) + radical(Allyl_S) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C][CH]C=C([O])[O](20902)',
    structure = SMILES('[CH2][C][CH]C=C([O])[O]'),
    E0 = (671.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,180,632.503,632.525,632.534,632.558],'cm^-1')),
        HinderedRotor(inertia=(0.0083209,'amu*angstrom^2'), symmetry=1, barrier=(2.36235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12273,'amu*angstrom^2'), symmetry=1, barrier=(112.985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107818,'amu*angstrom^2'), symmetry=1, barrier=(30.6185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17207,0.0601104,-6.12956e-05,3.15791e-08,-6.46474e-12,80838.5,26.2199], Tmin=(100,'K'), Tmax=(1181.72,'K')), NASAPolynomial(coeffs=[13.568,0.0181506,-8.03355e-06,1.53083e-09,-1.0775e-13,77908.8,-35.6531], Tmin=(1181.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(CCJ2_triplet) + radical(RCCJ) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][C]=[C]C[C]([O])[O](20903)',
    structure = SMILES('[CH2][C]=[C]C[C]([O])[O]'),
    E0 = (871.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360,370,350,180,1348.1,1350.73,1354.69],'cm^-1')),
        HinderedRotor(inertia=(0.179722,'amu*angstrom^2'), symmetry=1, barrier=(4.13216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180653,'amu*angstrom^2'), symmetry=1, barrier=(4.15356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182289,'amu*angstrom^2'), symmetry=1, barrier=(4.19117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37593,0.0716396,-0.000135847,1.40597e-07,-5.41183e-11,104866,30.1056], Tmin=(100,'K'), Tmax=(848.299,'K')), NASAPolynomial(coeffs=[1.09225,0.0397771,-2.08006e-05,4.04706e-09,-2.79415e-13,106109,38.4684], Tmin=(848.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(871.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ) + radical(Cds_S) + radical(CCOJ) + radical(Allyl_P) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C=CC=C([O])[O](20904)',
    structure = SMILES('[CH2]C=CC=C([O])[O]'),
    E0 = (92.9473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17348,0.0516073,-2.263e-05,-1.41036e-08,1.09533e-11,11290.3,23.9707], Tmin=(100,'K'), Tmax=(958.866,'K')), NASAPolynomial(coeffs=[15.9093,0.0152146,-4.93223e-06,8.69082e-10,-6.23405e-14,7311.42,-52.5143], Tmin=(958.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.9473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C([O])([O])[CH][C]=C(18606)',
    structure = SMILES('[CH2]C([O])([O])[CH][C]=C'),
    E0 = (610.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862103,0.0796233,-0.000133292,1.27933e-07,-4.7867e-11,73490.6,25.5549], Tmin=(100,'K'), Tmax=(820.903,'K')), NASAPolynomial(coeffs=[4.5812,0.0397071,-2.05311e-05,4.01678e-09,-2.80027e-13,73614.4,12.8191], Tmin=(820.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)(O)OJ) + radical(C=CCJCO) + radical(Cds_S) + radical(CC(C)(O)OJ) + radical(CJC(O)2C)"""),
)

species(
    label = 'C=C1[CH]CC1([O])[O](20905)',
    structure = SMILES('C=C1[CH]CC1([O])[O]'),
    E0 = (321.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.6691,-0.0185296,0.000118529,-1.12711e-07,2.5524e-11,38294.3,-22.9893], Tmin=(100,'K'), Tmax=(1761.06,'K')), NASAPolynomial(coeffs=[79.9043,0.0110857,-6.39334e-05,1.58861e-08,-1.18169e-12,-14364.4,-472.074], Tmin=(1761.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)(O)OJ) + radical(Allyl_S) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C=C=CC[C]([O])[O](18572)',
    structure = SMILES('C=C=CC[C]([O])[O]'),
    E0 = (420.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,540,610,2055,360,370,350,180,180,2398.26,2398.36],'cm^-1')),
        HinderedRotor(inertia=(0.193899,'amu*angstrom^2'), symmetry=1, barrier=(4.45812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57058,'amu*angstrom^2'), symmetry=1, barrier=(36.1108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50806,0.0642762,-0.000102925,1.0292e-07,-4.03372e-11,50673.6,27.4568], Tmin=(100,'K'), Tmax=(812.921,'K')), NASAPolynomial(coeffs=[2.15631,0.0403657,-2.05715e-05,4.02718e-09,-2.81671e-13,51252.9,28.6747], Tmin=(812.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][C][O](2319)',
    structure = SMILES('[O][C][O]'),
    E0 = (504.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([786.798,787.166,787.343],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.36312,0.000897068,5.83704e-06,-5.0699e-09,1.05455e-12,60693,5.58445], Tmin=(100,'K'), Tmax=(1870.77,'K')), NASAPolynomial(coeffs=[6.74008,0.00434757,-3.77128e-06,7.92205e-10,-5.4643e-14,58310.5,-11.3625], Tmin=(1870.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(CH2_triplet) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[CH2][C]=C[CH2](18777)',
    structure = SMILES('[CH2][C]=C[CH2]'),
    E0 = (510.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.997629,'amu*angstrom^2'), symmetry=1, barrier=(87.0906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0766176,'amu*angstrom^2'), symmetry=1, barrier=(26.8011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62092,0.0235622,5.21e-06,-2.26288e-08,1.0095e-11,61507.3,14.9566], Tmin=(100,'K'), Tmax=(985.296,'K')), NASAPolynomial(coeffs=[8.77433,0.0145496,-5.37933e-06,9.84612e-10,-6.99436e-14,59519.6,-18.5722], Tmin=(985.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC[C]([O])[O](20861)',
    structure = SMILES('[CH][C]=CC[C]([O])[O]'),
    E0 = (852.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,360,370,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43032,0.0692108,-0.000118391,1.23128e-07,-4.85611e-11,102622,30.265], Tmin=(100,'K'), Tmax=(838.443,'K')), NASAPolynomial(coeffs=[-0.222012,0.046687,-2.36968e-05,4.58127e-09,-3.16664e-13,103967,44.3189], Tmin=(838.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(852.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C[C]([O])[O](20853)',
    structure = SMILES('C#C[CH]C[C]([O])[O]'),
    E0 = (568.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,750,770,3400,2100,2175,525,3025,407.5,1350,352.5,186.946,994.087,994.979,999.961],'cm^-1')),
        HinderedRotor(inertia=(0.0642888,'amu*angstrom^2'), symmetry=1, barrier=(1.5289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0720289,'amu*angstrom^2'), symmetry=1, barrier=(1.70382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0766339,'amu*angstrom^2'), symmetry=1, barrier=(1.85403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38709,0.068768,-0.000121138,1.17164e-07,-4.23444e-11,68471.3,26.9609], Tmin=(100,'K'), Tmax=(884.473,'K')), NASAPolynomial(coeffs=[3.05765,0.0339367,-1.58081e-05,2.905e-09,-1.93158e-13,69242.7,25.1378], Tmin=(884.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(CCOJ) + radical(Sec_Propargyl) + radical(Cs_P)"""),
)

species(
    label = 'C=[C]C=C(18778)',
    structure = SMILES('C=[C]C=C'),
    E0 = (292.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.58595,'amu*angstrom^2'), symmetry=1, barrier=(36.4641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60774,0.0234239,8.49718e-06,-3.01392e-08,1.43429e-11,35286.5,12.5896], Tmin=(100,'K'), Tmax=(918.282,'K')), NASAPolynomial(coeffs=[9.56656,0.0118206,-3.11021e-06,4.74837e-10,-3.20206e-14,33219.7,-24.6845], Tmin=(918.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C[CH][C]([O])[O](18580)',
    structure = SMILES('C=[C]C[CH][C]([O])[O]'),
    E0 = (694.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,360,370,350,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,180,180,1556.16,1556.18,1556.23],'cm^-1')),
        HinderedRotor(inertia=(0.179002,'amu*angstrom^2'), symmetry=1, barrier=(4.1156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178957,'amu*angstrom^2'), symmetry=1, barrier=(4.11457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00239476,'amu*angstrom^2'), symmetry=1, barrier=(4.1156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55394,0.0641444,-0.000106618,1.08289e-07,-4.24313e-11,83551.7,31.647], Tmin=(100,'K'), Tmax=(824.246,'K')), NASAPolynomial(coeffs=[1.57624,0.0404595,-2.06092e-05,4.0203e-09,-2.80065e-13,84348.9,36.402], Tmin=(824.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cs_P) + radical(CCOJ) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH]C=CC[C]([O])[O](20855)',
    structure = SMILES('[CH]C=CC[C]([O])[O]'),
    E0 = (614.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66753,0.0408788,-1.82674e-05,2.74675e-09,-9.16583e-14,73912.5,21.5721], Tmin=(100,'K'), Tmax=(2766.07,'K')), NASAPolynomial(coeffs=[44.2147,-0.00506249,-2.37689e-07,6.03775e-11,1.19081e-15,46625.2,-224.076], Tmin=(2766.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_P) + radical(AllylJ2_triplet) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=[C]CC[C]([O])[O](18582)',
    structure = SMILES('[CH]=[C]CC[C]([O])[O]'),
    E0 = (741.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,360,370,350,180,180,1964.91,1965.18],'cm^-1')),
        HinderedRotor(inertia=(0.12083,'amu*angstrom^2'), symmetry=1, barrier=(2.77812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121282,'amu*angstrom^2'), symmetry=1, barrier=(2.78852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120963,'amu*angstrom^2'), symmetry=1, barrier=(2.78117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4214.4,'J/mol'), sigma=(6.97559,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.28 K, Pc=28.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28187,0.0729175,-0.000131969,1.35104e-07,-5.20697e-11,89235,30.5305], Tmin=(100,'K'), Tmax=(841.972,'K')), NASAPolynomial(coeffs=[1.35808,0.041739,-2.15235e-05,4.18441e-09,-2.89526e-13,90314.5,36.6625], Tmin=(841.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P) + radical(CCOJ) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = '[CH][C]=CCC([O])[O](20856)',
    structure = SMILES('[CH][C]=CCC([O])[O]'),
    E0 = (647.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10934,0.0394501,-1.68094e-05,2.13179e-09,-1.14362e-14,77808.1,19.712], Tmin=(100,'K'), Tmax=(2796.54,'K')), NASAPolynomial(coeffs=[50.6318,-0.011229,1.86462e-06,-2.91431e-10,2.38529e-14,45584.3,-263.666], Tmin=(2796.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC[C]([O])O(20859)',
    structure = SMILES('[CH][C]=CC[C]([O])O'),
    E0 = (626.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10863,0.0714467,-0.000104151,9.52798e-08,-3.51671e-11,75491.8,31.2545], Tmin=(100,'K'), Tmax=(816.879,'K')), NASAPolynomial(coeffs=[4.36771,0.0404497,-1.96185e-05,3.75559e-09,-2.59692e-13,75461.1,19.2616], Tmin=(816.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_P) + radical(Cds_S) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]CC=C([O])[O](18571)',
    structure = SMILES('C=[C]CC=C([O])[O]'),
    E0 = (243.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,350,440,435,1725,418.018,418.045,418.048,418.071],'cm^-1')),
        HinderedRotor(inertia=(0.107796,'amu*angstrom^2'), symmetry=1, barrier=(13.371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107821,'amu*angstrom^2'), symmetry=1, barrier=(13.3716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27347,0.0581383,-5.45988e-05,2.65565e-08,-5.21137e-12,29396.5,25.5517], Tmin=(100,'K'), Tmax=(1218.9,'K')), NASAPolynomial(coeffs=[12.5083,0.0212695,-9.22758e-06,1.74115e-09,-1.21688e-13,26657.6,-30.8741], Tmin=(1218.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=C=C[CH]C([O])[O](20906)',
    structure = SMILES('C=C=C[CH]C([O])[O]'),
    E0 = (332.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29276,0.0670164,-9.85606e-05,9.25329e-08,-3.56207e-11,40059.7,24.7416], Tmin=(100,'K'), Tmax=(775.397,'K')), NASAPolynomial(coeffs=[4.26391,0.0390447,-1.99887e-05,3.94787e-09,-2.78785e-13,39979.1,13.6145], Tmin=(775.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C1CC1([O])[O](20907)',
    structure = SMILES('C=[C]C1CC1([O])[O]'),
    E0 = (431.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73134,0.0436482,-2.38843e-05,4.67763e-09,2.12414e-14,51985.1,24.3013], Tmin=(100,'K'), Tmax=(1428.35,'K')), NASAPolynomial(coeffs=[11.8434,0.0228355,-9.90945e-06,1.83383e-09,-1.24912e-13,48330.8,-30.7688], Tmin=(1428.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)(O)OJ) + radical(Cds_S) + radical(CC(C)(O)OJ)"""),
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
    E0 = (633.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (633.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (787.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (633.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (762.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (751.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (837.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (807.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (826.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (792.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (881.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (793.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (764.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1026.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (883.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1083.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (696.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (878.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (641.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (633.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1050.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1064.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (795.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (797.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (852.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (828.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (886.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (680.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (762.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (656.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (656.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (636.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=C([O])[O](1172)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH2][C]=CC=C([O])[O](20892)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(90.798,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 85.8 to 90.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH2]C#CC[C]([O])[O](20893)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])[O](1172)', '[CH][C]=C(18825)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.123469,'m^3/(mol*s)'), n=2.00579, Ea=(59.5966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 54.2 to 59.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C]([O])[O](2304)', 'C#C[CH2](17441)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH2][C]=C[CH]C([O])[O](20894)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.27137e+08,'s^-1'), n=1.53496, Ea=(117.681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=[C]C[C]([O])[O](20895)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=[C]CC([O])[O](20896)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH2][CH][CH]C=C([O])[O](20897)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH2][C]=C[CH][C]([O])O(20898)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.64e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/OneDe] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[C]=[C]C[C]([O])[O](20899)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=[C]C[C]([O])O(20900)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C[C][CH]C=C([O])[O](20901)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]([O])[O](2304)', '[CH][C]=C(18825)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH2][C][CH]C=C([O])[O](20902)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH2][C]=[C]C[C]([O])[O](20903)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH2]C=CC=C([O])[O](20904)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH2]C([O])([O])[CH][C]=C(18606)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=C1[CH]CC1([O])[O](20905)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=C=CC[C]([O])[O](18572)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C][O](2319)', '[CH2][C]=C[CH2](18777)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[CH][C]=CC[C]([O])[O](20861)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(8)', 'C#C[CH]C[C]([O])[O](20853)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C][O](2319)', 'C=[C]C=C(18778)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C[CH][C]([O])[O](18580)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH]C=CC[C]([O])[O](20855)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC[C]([O])[O](18582)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH][C]=CCC([O])[O](20856)'],
    products = ['[CH2][C]=CC[C]([O])[O](18581)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['[CH][C]=CC[C]([O])O(20859)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.854,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6Hall;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=[C]CC=C([O])[O](18571)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=C=C[CH]C([O])[O](20906)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]=CC[C]([O])[O](18581)'],
    products = ['C=[C]C1CC1([O])[O](20907)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

network(
    label = '4005',
    isomers = [
        '[CH2][C]=CC[C]([O])[O](18581)',
    ],
    reactants = [
        ('C=C([O])[O](1172)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4005',
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


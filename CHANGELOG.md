
version [2.4.2]
  - fixes bug for time variable light curve

version [2.4.1]
  - fixes spectrum caching
    (fixes segfault when simulating more than ten sources)

version [2.4.0]
  - setting ARF to 0.0 for any bit NOT >=0
    (fixes SIXTE crash if an ARF contained NULL values)
  - increases normalization criterion in loadRMF to 10%

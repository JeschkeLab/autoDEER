# Logging: An EPR implementation
When building any large piece of software one of the most important and first features that need to be created in a logging function. This program is no different, and will also have one added. Since this software is primarily written in python we propose using the integrated and built logging module.

The structure of this logging will roughly follow that of the module in general, with a separate hardware log to the core modules.

## Logging Levels
One of the main features of the python logging modules is the accessibility of multiple logging levels. There are five default levels that we will aim to use however we might add some extra at a later date.
1) DEBUG: This is the lowest level and generally shouldn't be printed to the log file unless there is an issue that needs to be debuged.
2) INFO: This will be the lowest level that is always recorded. Here we note down everything that happens: An experiment starts, a parameter changes etc… This will could also be used by the user to retrospectively investigate what has happened.
3) WARNING: This is xthe lowest level that suggests a problem. Here we are recording something that a problem might occur. This could be hardware related or it could be that there is a high degree of uncertainty in a calculated parameter. These warning should also be passed to the user in a more accessible format, maybe a pop-up or so. One of the most likely causes of this is that the temperature has risen outside of an acceptable region. The experiment will be saved but not stopped.
4) ERROR: Something has gone wrong. This **MUST** be passed on to the user and all hardware should be sent to a safe state. I.e attenuators set to to maximum, AWG stoped etc… The experiment will be saved and stopped.
5) CRITICAL: Something serious has gone wrong and the safety of the equipment is danger. This could be that defence pulses are missing. That too much power was detected. ALL Attenuators will immediately be set to maximum **and** the amplifier will be switched to stand by or off. This will be reported to the user at the desktop **and** a alert will be sent to either the user or the administrator by the selected communication channel (email or slack). After a critical alert all hardware should be carefully checked before anything is restarted. The experiment will be stopped, only the last autosave will be kept

## Logs
There will be two logs. A core _log_ and a _hardware_ log. These will both be saved in the appropriate tmp folder. On unix systems this /tmp/autoDEER. All logs will be encoded in "UTF-8" per the standards for this program with exception of file names that will be encoded with "ASCII". Since experiments and software are often run for longer than a day even weeks, both the date and time will be recorded for each log.

Logs are created with these two commands, in the main script that runs the whole core.:
```python
import logging

logging.basicConfig(filename="core.log",format='%(levelname)s: %(asctime)s %(message)s', encoding="utf-8",level.logging.DEBUG)
logging.basicConfig(filename="hardware.log",format='%(levelname)s: %(asctime)s %(message)s', encoding="utf-8",level.logging.DEBUG)

```

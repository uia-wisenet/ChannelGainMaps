function ch_expNumber = whichExp()
str_stack = dbstack;
theClass = str_stack(2).file(1:end-1);
ch_expNumber = char(extractAfter(...
    str_stack(2).name, theClass + "experiment_"));


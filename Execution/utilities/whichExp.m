function ch_expNumber = whichExp()
ch_expNumber = '';
try
    obj = evalin('caller', 'obj');
    str_stack = dbstack;
    theClass = str_stack(2).file(1:end-1);
    ch_expNumber = char(extractAfter(...
    str_stack(2).name, strcat(theClass, obj.funbasename)));
catch me
    switch me.identifier
        case 'MATLAB:UndefinedFunction'
            warning 'whichExp must be called within a method of a class inheriting from ExperimentFunctionSet'
        case 'MATLAB:badsubscript'
            warning 'whichExp must be called within a method of a class inheriting from ExperimentFunctionSet'
    end
end

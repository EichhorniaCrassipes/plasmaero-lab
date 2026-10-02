function hantek(id,waveform,params)
    arguments
        id {mustBeText}
        waveform {mustBeMember(waveform, {'DC'})}
        params.channel (1,:) double
        params.amplitude (1,:) double
        params.frequency (1,:) double
    end
    persistent fg;
    try
        if isempty(fg)
            fg = visadev(id);
            fprintf('Connected to Hantek HDG3000B.\n');

            idn = writeread(fg, '*IDN?');
            fprintf('Device Identity: %s\n', strtrim(idn));
        end
    catch ME
        error('Communication failed: %s', ME.message);
    end
    try
        switch waveform
            case 'DC'
                if ~isfield(params, 'channel') || ~isfield(params, 'voltage')
                    error('Missing parameters params.channel and params.voltage')
                end
                ch = params.channel;
                v_dc = params.voltage;
                % Step A: Set the function wave shape to DC
                write(fg, sprintf('SOURce%d:FUNCtion DC', ch));

                % Step B: Apply the DC offset voltage level
                write(fg, sprintf('SOURce%d:VOLTage:OFFSet %.3f', ch, v_dc));

                % Step C: Enable the Output Channel
                write(fg, sprintf('OUTPut%d ON', ch));

                fprintf('Channel %d is now active and outputting %.2f V DC.\n', ch, v_dc);

            otherwise
                error('No mode "%s" in function. Maybe you mean "DC"?', waveform)
        end
    catch ME
        error('Communication failed: %s', ME.message);
    end
end

%% function usage example
visaAddress = 'USB0::0x0483::0x5740::CN2225055000657::0::INSTR';

dc_params.channel = 1;
dc_params.voltage = 2.5;

hantek(visaAddress, 'DC', dc_params);


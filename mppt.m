%% 
% *Modified P&O Maximum Power Point Tracking*
% Copyright 2022 Henri Salonen
% Solar cell measurements designed jointly with Arto Hiltunen
%This program is free software: you can redistribute it and/or
%modify it under the terms of the GNU General Public License 
%as published by the Free Software Foundation, either version 3
%of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.

%You should have received a copy of the 
%GNU General Public License along with this program. 
%If not, see <https://www.gnu.org/licenses/>.

% Designed as a lightweight MPP tracking tool for 
% experimental solar cell research.
% Requirements: NI-VISA, MATLAB 2006a onwards, compliant SMU
% Tested on MATLAB 2021a, Keithley 2450 SMU

function mppt(name, u_start, u_stop, a_sample)
%% 
% *Function call parameters:*

%       name        : savename for data, ie. 'example_name'
%       u_start     : initial voltage of sweep 1 (V)
%       u_stop      : final voltage of sweep 1 (V)
%       a_sample    : area of sample (mm^2)
%%
% *User settings*

meas_type = 3;          % Type of measurement, options:
                        % 1 : measure only both IV-sweeps
                        % 2 : measure MPPT with a single IV-sweep
                        % 3 : measure MPPT with both IV-sweeps
%           [H, M, S]                                                           
meas_time = [0, 1, 0];  % Maximum length of measurement,                    

% SMU parameters
current_range = ' 100E-3'; % Options:
                           % ':AUTO ON'    - Autorange
                           % ' 100E-6'     - 100uA
                           % ' 10E-3'      - 10mA
                           % ' 100E-3'     - 100mA
                           % ' 1'          - 1A

voltage_range = ' 2';      % Options:
                           % ':AUTO ON'    - Autorange
                           % ' 2'          - 2 V
                           % ' 20'         - 20 V
                                
current_limit = ' 100E-3'; % Options, must be equal to
                           % or higher than current range:
                           % ' 10E-3'      - 10mA
                           % ' 100E-3'     - 100mA
                           % ' 1'          - 1A

visa_addr = 'USB0::0x05E6::0x2450::04386242::INSTR'; % SMU address
                                
% IV sweep parameters
IV_step = 0.01;         % V / step
IV_rate = 0.2;          % Rate
pre_sweep_pause = 3;    % Voltage hold / s

% MPPT parameters
phase1_step = 0.1;      % Initial step size    / V
phase2_step = 0.01;     % Step size after MPP  / V
v_interval = 2;         % Tdwell length        / s
v_start = 0;            % Initial voltage      / V
p_lamp = 0.1;           % Power of lamp        / W/cm^2
m_start_dir_up = true;  % true - start near 0 
                        % false - start near Voc
v_from_sweep = true;    % true - v_start from first IV


%%
% *Script Starts*

% Initial verification of name variable
if(~isvarname(name))
    disp(['Name invalid! Name cannot start with'...
        ' a number or contain certain special characters'])
    return
end

% Set up parameters for the IV-sweep
% Calculate the number of steps
num_steps = abs(u_stop-u_start) ./ IV_step +1;
param = [num2str(u_start) ', ' num2str(u_stop)...
    ', ' num2str(num_steps) ', ' num2str(IV_rate)];
param_rev = [num2str(u_stop) ', ' num2str(u_start)...
    ', ' num2str(num_steps) ', ' num2str(IV_rate)];

% Set up names for data output and files
name_iv = [name '_sweep'];
name_mppt = [name '_mppt'];
filename_mppt = [name_mppt '.csv'];

% Set up variables for data and control
m_data_t = [0,0];
m_data_eff = [0,0];     
% Measurement direction TRUE: Vn > Vn-1 , FALSE: Vn < Vn-1
m_dir_up = m_start_dir_up; 
phase_one = true;
op_t = 0;
% Default value for power to be compared against
p_previous = -10;
v_now = v_start;
meas_time = sum(meas_time.*[3600, 60, 1]);
smu_open = false;
count = floor(25*v_interval);
half_count = ceil(count/2);

%% 
% *SMU connection and settings*
instrreset

try % Connect to the SMU
    visaSMU = visa('ni',visa_addr);
    visaSMU.InputBufferSize = 100000;
    visaSMU.Timeout = 15;
    visaSMU.ByteOrder = 'littleEndian';
    fopen(visaSMU);
    smu_open = true;
catch
    error( 'SMU init failed' );
end
% Reset settings on the SMU
fprintf(visaSMU, '*RST');                           
% Set to act as voltage source
fprintf(visaSMU, 'SOUR:FUNC VOLT');                 
% Set voltage delay to 0
fprintf(visaSMU, 'SOUR:VOLT:DEL 0');                
% Automatic zeroing once before measurement
fprintf(visaSMU, 'SENS:AZER:ONCE');            
% Measure DC current
fprintf(visaSMU, 'SENS:FUNC "CURR"');               
% Set 4-wire measurement
fprintf(visaSMU, 'SENS:CURR:RSEN ON');              
% Set the number of power line cycles
fprintf(visaSMU, 'SENS:NPLC 1');                   
fprintf(visaSMU, ['SENS:CURR:RANG' current_range]); % Set I range
fprintf(visaSMU, ['SOUR:VOLT:RANG' voltage_range]); % Set V range
fprintf(visaSMU, ['SOUR:VOLT:ILIM' current_limit]); % Set I limit

%% 
% *IV sweep 1*
fprintf(visaSMU, 'OUTP ON');
fprintf(visaSMU, ['SOUR:VOLT ', num2str(u_start,'%2.3f')]);
% Hold sweep starting voltage in preparation
pause(pre_sweep_pause);            
% Set sweep function measurement
fprintf(visaSMU, ['SOUR:SWE:VOLT:LIN ' param]);
% Start measurement
fprintf(visaSMU, ':INIT');
% SMU: Wait for completion, disregard commands until done
fprintf(visaSMU, '*WAI');
% Pause matlab, wait until sweep is finished
pause(num_steps*IV_rate -1);                    

% Read sweep data, separate and process
data_s1 = query(visaSMU, ['TRAC:DATA? 1, '...
    num2str(num_steps) ', "defbuffer1",  SOUR, READ']);
data_s1 = split(data_s1, ',');
data_s1 = cellfun(@str2num,data_s1);
data_v_s1 = data_s1(1:2:end);
data_i_s1 = data_s1(2:2:end);
data_iv_s1 = [data_v_s1 data_i_s1];
assignin('base',[name_iv '1'], data_iv_s1);
csvwrite([name_iv '1.csv'], data_iv_s1);
fig_iv = figure(1);
clf(fig_iv);
% Draw sweep 1 IV-data
figure(1);
plot(data_iv_s1(:,1), data_iv_s1(:,2));

fprintf(visaSMU, 'TRAC:CLEAR "defbuffer1"');

%%
% *IV sweep 2*

if( meas_type == 1 || meas_type == 3)
    fprintf(visaSMU, 'OUTP ON');
    fprintf(visaSMU, ['SOUR:VOLT ', num2str(u_stop,'%2.3f')]);
    pause(pre_sweep_pause);                             
    % Measure and wait
    fprintf(visaSMU, ['SOUR:SWE:VOLT:LIN ' param_rev]); 
    fprintf(visaSMU, ':INIT');                          
    fprintf(visaSMU, '*WAI');                           
    pause(num_steps*IV_rate -1);                           
    
    % Read sweep data, separate and process
    data_s2 = query(visaSMU, ['TRAC:DATA? 1, '...
        num2str(num_steps) ', "defbuffer1",  SOUR, READ']);
    data_s2 = split(data_s2, ',');
    data_s2 = cellfun(@str2num,data_s2);
    data_v_s2 = data_s2(1:2:end);
    data_i_s2 = data_s2(2:2:end);
    data_iv_s2 = [data_v_s2 data_i_s2];
    assignin('base',[name_iv '2'], data_iv_s2);
    csvwrite([name_iv '2.csv'], data_iv_s2);
    % Draw sweep 2 IV-data
    figure(1);
    hold on;
    plot(data_iv_s2(:,1), data_iv_s2(:,2));
    hold off;
    fprintf(visaSMU, 'TRAC:CLEAR "defbuffer1"');
end

%% 
% *Figures of merit*

% Sweep 1
disp('1st Sweep FoM');
% Find all negative voltage values
Isc_s1 = find(data_v_s1 <= 0);      
% The closest to zero is Isc
[~, ind_i] = max(data_v_s1(Isc_s1));            
Isc_s1 = data_i_s1(Isc_s1(ind_i));
% Find all negative current values
Voc_s1 = find(data_i_s1 <= 0);      
% The closest to zero is Voc
[~, ind_v] = max(data_i_s1(Voc_s1));            
Voc_s1 = data_v_s1(Voc_s1(ind_v));
% Pmax is the minimum, as current is negative
[Pmax_s1, ind_v_mpp] = min(data_i_s1 .* data_v_s1);      
FF_s1 = Pmax_s1 / (Voc_s1 *Isc_s1) * 100;
eff_s1 = Pmax_s1 / (p_lamp * (a_sample/100)) *(-1)*100;
Isc_s1 = (Isc_s1/a_sample)*100*1000*(-1);

columntitle = {'efficiency (%) ', 'FF ',...
    'Isc (mA/cm^2) ', 'Voc (V) ' ,'area (mm^2) ' };
result = [eff_s1 FF_s1 Isc_s1 Voc_s1 a_sample];
FoM = [columntitle ; num2cell(result)];
disp(FoM)
assignin('base', [name '_FoM_s1'], FoM)

% Sweep 2
if( meas_type == 1 || meas_type == 3 )
    disp('2nd Sweep FoM');
    Isc_s2 = find(data_v_s2 <= 0);              
    [~, ind_i] = max(data_v_s2(Isc_s2));           
    Isc_s2 = data_i_s2(Isc_s2(ind_i));
    Voc_s2 = find(data_i_s2 <= 0);              
    [~, ind_v] = max(data_i_s2(Voc_s2));            
    Voc_s2 = data_v_s2(Voc_s2(ind_v));
    [Pmax_s2, ind_v_mpp] = min(data_i_s2 .* data_v_s2);      
    FF_s2 = Pmax_s2 / (Voc_s2 *Isc_s2) * 100;
    eff_s2 = Pmax_s2 / (p_lamp * (a_sample/100)) *(-1)*100;
    Isc_s2 = (Isc_s2/a_sample)*100*1000*(-1);
    
    result2 = [eff_s2 FF_s2 Isc_s2 Voc_s2 a_sample];
    FoM_s2 = [columntitle ; num2cell(result2)];
    disp(FoM_s2)
    assignin('base', [name '_FoM_s2'], FoM_s2)
end

%% 
% *MPP tracking*
if( meas_type == 2 || meas_type == 3 )
% Require confirmation of a working reading
meas_cont = questdlg('Continue with measurement?',...
    'Check IV curve', 'Yes', 'No', 'Yes');

meas_cont = strcmp(meas_cont, 'Yes');
end
if( (meas_type == 2 || meas_type == 3) && meas_cont)
        
    % Ensure SMU settings
    fprintf(visaSMU, 'SOUR:FUNC VOLT'); 
    fprintf(visaSMU, 'SOUR:VOLT:DEL 0');
    fprintf(visaSMU, 'SENS:AZER:ONCE'); 
    fprintf(visaSMU, 'SENS:FUNC "CURR"'); 
    fprintf(visaSMU, 'SENS:CURR:RSEN ON'); 
    fprintf(visaSMU, 'SENS:NPLC 1'); 
    fprintf(visaSMU, 'TRAC:CLEAR "defbuffer1"');
    % Calculates the number of datapoints (NPLC 1, read back on)
    fprintf(visaSMU, ['COUNT ' num2str(count)]); 

    % Create the figure and controls for monitoring
    fig_m = uifigure;
    fig_m.Position = [100 100 650 400];
    fig_m.Name = [name ' MPPT'];
    axis_m = uiaxes(fig_m);
    axis_m.Position = [10 75 500 300];
    axis_m.XGrid = 'on';
    axis_m.YGrid = 'on';
    cb_m = uicheckbox(fig_m, 'Text',...
        'Stop measurement', 'Value', 0);
    cb_m.Position = [10 10 160 22];
    plot_m = plot(axis_m, m_data_t, m_data_eff);
    plot_m.XDataSource = 'm_data_t';
    plot_m.YDataSource = 'm_data_eff';
    c_s_index = '1';
    data_t_now = [];
    data_v_now = [];
    data_i_now = [];

    % If v_from_sweep is true, set v_start as v_mpp
    % and start with small steps
    if( v_from_sweep )
        phase_one = false;
        v_start = data_v_s1(ind_v_mpp);
        v_now = v_start;
    end
    
    % Apply the defined starting voltage
    fprintf(visaSMU, 'OUTP ON');
    fprintf(visaSMU, ['SOUR:VOLT ', num2str(v_start,'%2.3f')]);

    start_time = tic;
    elapsed_time = toc(start_time);

    % Main loop
    while(cb_m.Value == 0 && elapsed_time<meas_time)
    % Measure v_interval*25 (NPLC 1 & read back on)
        fprintf(visaSMU, 'TRAC:TRIG "defbuffer1"');
        pause(0.005)

        op_start_t = tic;
        % Prepare for communication errors
        try
            c_e_index = strtrim(query(visaSMU,...
                'TRAC:ACT:END? "defbuffer1"'));
        catch E_q
            e_log = E_q;
            assignin('base','e_log', e_log);
            assignin('base','backup_data_t', data_t_now);
            assignin('base','backup_data_v', data_v_now);
            assignin('base','backup_data_i', data_i_now);
            disp('Query error');
            smu_open = false;
            instrreset;
            break;
        end
        % Gather data and parse
        try
            data_now = query(visaSMU, ['TRAC:DATA? ' c_s_index...
                ',' c_e_index ',"defbuffer1", SOUR, READ, TST']);
        catch E_q
            e_log = E_q;
            assignin('base','e_log', e_log);
            assignin('base','backup_data_t', data_t_now);
            assignin('base','backup_data_v', data_v_now);
            assignin('base','backup_data_i', data_i_now);
            disp('Query error');
            instrreset;
            smu_open = false;
            break;
        end    
        c_s_index = num2str(floor(str2double(c_e_index)+1));

        % Every 1h, clear buffer
        if( str2double(c_s_index) > 3600) 
            fprintf(visaSMU, 'TRAC:CLEAR "defbuffer1"');
            c_s_index = '1';
            dlmwrite('mppt_backup.csv', [m_data_eff data_v_now...
                data_i_now m_data_t],...
                'delimiter',',','precision',12);
        end
        data_now = split(data_now, ',');
        data_t_now = [data_t_now ; data_now(3:3:end)];
        data_now = cellfun(@str2double, data_now);
        assignin('base', 'data_all', data_now);
        data_v_now = [data_v_now ; data_now(1:3:end)];
        data_i_now = [data_i_now ; data_now(2:3:end)];

        % Calculate P and PCE, generate rough timestamps
        p_now = data_v_now.*data_i_now;
        m_data_eff = (p_now ./ (p_lamp*(a_sample/100))) ;
        m_data_eff = m_data_eff .* (-1);
        m_data_t = (1:length(m_data_eff)).*0.04;

        % Check direction
        % Average the latest samples
        p_current = (-1)*sum(p_now(end-(half_count)...
            +1:end))/half_count;

        % If Pn < Pn-1, change direction
        if(p_current < p_previous)
            if(phase_one) % Switch to phase two
                phase_one = false;
            end
            m_dir_up = ~m_dir_up;
        end
        
        % Update Pn-1 as Pn for next iteration
        p_previous = p_current;
        
        % If in phase one, larger steps
        if(phase_one)
            if( m_dir_up)
                v_now = v_now + phase1_step;
            else
                v_now = v_now - phase1_step;
            end
        else
            % Else use smaller steps
            if(m_dir_up)
                v_now = v_now + phase2_step;
            else
                v_now = v_now - phase2_step;
            end
            % Periodic flush to prevent buffer overflow
            flushinput(visaSMU);
            flushoutput(visaSMU);
        end

        fprintf(visaSMU, ['SOUR:VOLT ',...
            num2str(v_now,'%2.3f')]);

        % Refresh the data for the plot
        assignin('base', 'm_data_t', m_data_t);
        assignin('base', 'm_data_eff', m_data_eff);
        assignin('base', 'm_data_i', data_i_now);
        assignin('base', 'm_data_v', data_v_now);
        
        refreshdata(plot_m);
        drawnow();
        elapsed_time = toc(start_time);
        op_t = toc(op_start_t);
        assignin('base', 'op_t', op_t);
    end
    
    if(smu_open)
    fprintf(visaSMU, 'OUTP OFF');
    fclose(visaSMU);
    smu_open = false;
    end
    
    disp('Parsing timestamps, please wait');
    % Parse actual timestamps
    % Format of timestamp
    f_tst = 'MM/dd/yyyy HH:mm:ss.SSSSSSSSS'; 
    tst = datetime(data_t_now, 'InputFormat',f_tst);
    tst = seconds(tst-tst(1));
    
    data_all = [m_data_eff data_v_now data_i_now tst];
    assignin('base', [name_mppt '_data'], data_all);
    disp('Writing csv data, please wait');
    dlmwrite(filename_mppt, data_all,...
        'delimiter',',','precision',12);

    msgbox('Measurement complete', 'Done', 'help');
        
end

%%
% *Shutdown procedures*
if(smu_open)
    fprintf(visaSMU, 'OUTP OFF');
    fclose(visaSMU);
end
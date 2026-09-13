classdef RequestHandler < handle
% REQUESTHANDLER  Minimal MATLAB client for the BOPTEST REST API.
%
% Usage:
%   handler = RequestHandler('http://127.0.0.1', 80);
%   tid     = handler.deploy_test('multizone_residential_hydronic');
%   handler.set_scenario(start_time, warmup_period);
%   handler.set_step(900);
%   res     = handler.advance(overrides_struct);
%   kpi     = handler.get_kpi();
%   handler.stop_test(tid);

    properties
        IP                    % base URL, e.g. 'http://127.0.0.1:80'
        TestList    = {}      % ids of every test deployed by this handler
        CurrentTest = []      % id the other methods act on
    end

    methods
        function obj = RequestHandler(ip_address, port)
            if nargin < 1, ip_address = 'http://127.0.0.1'; end
            if nargin < 2, port = 80; end
            obj.IP = sprintf('%s:%d', ip_address, port);
        end

        function test_id = deploy_test(obj, case_name)
            url  = obj.IP + "/testcases/" + case_name + "/select";
            resp = webwrite(url, struct(), weboptions('Timeout', 60));
            test_id         = resp.testid;
            obj.TestList{end+1} = test_id;
            obj.CurrentTest = test_id;
            fprintf('Deployed: %s  (ID: %s)\n', case_name, test_id);
        end

        function set_step(obj, step_seconds)
            webwrite(obj.IP + "/step/" + obj.CurrentTest, ...
                     struct('step', step_seconds), ...
                     weboptions('RequestMethod', 'put'));
            fprintf('Simulation step: %d s.\n', step_seconds);
        end

        function response = set_scenario(obj, start_time, warmup_period, ...
                                         temperature_uncertainty, solar_uncertainty)
        % SET_SCENARIO  Set start time, warm-up and optional forecast noise.
        %
        %   BOPTEST splits these over two endpoints, and confusing them fails
        %   silently. /scenario carries the forecast uncertainty and knows
        %   nothing about start_time or warmup_period: it returns 200 and
        %   ignores them, leaving the case where deploy_test put it, at t = 0.
        %   The scenario window belongs to /initialize. Both are called here in
        %   the order BOPTEST's own client uses — scenario first, then
        %   initialize — and the resulting time is read back and checked.
        %
        %   temperature_uncertainty / solar_uncertainty accept 'low' | 'medium' |
        %   'high' and activate BOPTEST's forecast noise on dry-bulb temperature
        %   and global horizontal irradiation. An empty value is omitted from the
        %   request body, so that channel stays deterministic.

            if nargin < 4 || isempty(temperature_uncertainty), temperature_uncertainty = ''; end
            if nargin < 5 || isempty(solar_uncertainty),       solar_uncertainty       = ''; end

            opts = weboptions('RequestMethod', 'put', 'Timeout', 120);

            unc_tag = '';   % only reported when a noise channel is active
            if ~isempty(temperature_uncertainty) || ~isempty(solar_uncertainty)
                tu = temperature_uncertainty; if isempty(tu), tu = 'det'; end
                su = solar_uncertainty;       if isempty(su), su = 'det'; end
                unc_tag = sprintf(', T_unc=%s, S_unc=%s', tu, su);

                body = struct();
                if ~isempty(temperature_uncertainty)
                    body.temperature_uncertainty = temperature_uncertainty;
                end
                if ~isempty(solar_uncertainty)
                    body.solar_uncertainty = solar_uncertainty;
                end
                webwrite(obj.IP + "/scenario/" + obj.CurrentTest, body, opts);
            end

            fprintf('Setting scenario (t0=%d, warmup=%d s%s)... ', ...
                    start_time, warmup_period, unc_tag);
            response = webwrite(obj.IP + "/initialize/" + obj.CurrentTest, ...
                                struct('start_time',    start_time, ...
                                       'warmup_period', warmup_period), opts);

            % Read the time back. A window that never takes effect is invisible
            % otherwise, and every episode then silently runs on the same week.
            t0 = obj.response_time(response);
            if isnan(t0)
                fprintf('Done (start time not reported).\n');
            else
                fprintf('Done, initialised at t=%d s (day %.2f).\n', round(t0), t0/86400);
                if abs(t0 - start_time) > 1
                    warning('RequestHandler:startTime', ...
                            ['BOPTEST initialised at %g s, not the requested %g s. ' ...
                             'The episode would run on the wrong week.'], t0, start_time);
                end
            end
        end

        function t = response_time(~, response)
        % RESPONSE_TIME  Pull 'time' out of a REST response, NaN if absent.
            t = NaN;
            if ~isstruct(response), return; end
            if isfield(response, 'payload') && isstruct(response.payload)
                p = response.payload;
            else
                p = response;
            end
            if isfield(p, 'time') && isnumeric(p.time) && isscalar(p.time)
                t = double(p.time);
            end
        end

        function measurements = advance(obj, overrides)
            if nargin < 2, overrides = struct(); end
            resp         = webwrite(obj.IP + "/advance/" + obj.CurrentTest, ...
                                    overrides, weboptions('Timeout', 60));
            measurements = resp.payload;
        end

        function forecast = get_forecast(obj, point_names, horizon, interval)
        % GET_FORECAST  Fetch forecasts over a horizon [s] at a given interval [s].
            body = jsonencode(struct('point_names', {point_names}, ...
                                    'horizon', horizon, 'interval', interval));
            opts = weboptions('RequestMethod', 'put', ...
                              'MediaType', 'application/json', ...
                              'CharacterEncoding', 'UTF-8', 'Timeout', 60);
            resp = webwrite(obj.IP + "/forecast/" + obj.CurrentTest, body, opts);
            if isfield(resp, 'payload'), forecast = resp.payload; else, forecast = resp; end
        end

        function data = get_results(obj, point_names, start_time, final_time)
            body = struct('point_names', {point_names}, ...
                          'start_time',  start_time, 'final_time', final_time);
            opts = weboptions('RequestMethod', 'put', ...
                              'MediaType', 'application/json', 'Timeout', 120);
            fprintf('Fetching results [%d, %d]... ', start_time, final_time);
            data = webwrite(obj.IP + "/results/" + obj.CurrentTest, body, opts);
            fprintf('Done.\n');
        end

        function kpi = get_kpi(obj)
        % GET_KPI  Official BOPTEST KPIs, cumulative from start_time (warm-up excluded).
        %   ener_tot [kWh/m2], cost_tot [$/m2], emis_tot [kgCO2/m2],
        %   tdis_tot [K*h/zone], idis_tot [ppm*h/zone], time_rat [-].
            resp = webread(obj.IP + "/kpi/" + obj.CurrentTest, ...
                           weboptions('Timeout', 30));
            if isfield(resp, 'payload')
                kpi = resp.payload;
            else
                kpi = resp;
            end
        end

        function stop_test(obj, test_id)
            webwrite(obj.IP + "/stop/" + test_id, struct(), ...
                     weboptions('RequestMethod', 'put'));
            obj.TestList(strcmp(obj.TestList, test_id)) = [];
        end
    end
end

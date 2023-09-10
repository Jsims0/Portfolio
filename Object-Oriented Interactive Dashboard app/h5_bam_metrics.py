import json
from os.path import basename, join
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd
import dash
from XXXXXXXXX.app import app
from XXXXXXXXX.libs.module import TelescopeModule
from XXXXXXXXX.libs.custom_components import CardCol
from XXXXXXXXX.libs.data.H5_bam_provider import H5LocalProvider
from XXXXXXXXX.libs.utils.common import get_triggered_component
from XXXXXXXXX.libs.utils.url_b64 import encode_to_b64, decode_from_b64, hash_strlist
from XXXXXXXXX.libs.spatial_error_helper import get_jobs_list_table
from XXXXXXXXX.libs.utils.path import Platform
from XXXXXXXXX.libs.utils.qsub import QSubManager



def local_id(x): return f'H5-bam-metrics-{x}'
# dbc toggle is a radio group and the off value is an empty list and on value is just [True]
# we define these as constants here to make it more readable when we assign these values
TOGGLE_OFF = []
TOGGLE_ON = [True]


def create_layout(app, local_config, url):
    # local_config = app.XXXXXXXXX_config['h5_bam_metrics']
    datasource = app.XXXXXXXXX_databroker.get(H5LocalProvider)
    #datasource_bm = app.XXXXXXXXX_databroker.get(H5BasemountProvider)
    default_metric_values = local_config['default_metric_set']
    sample_size = local_config['sample_size']
    # chart type / xvar / y var / colorby, sizeby, animateby, groupby
    default_chart_control_values = ['scatter', 6, None, None, None, None, None]

    if app.XXXXXXXXX_download_enabled:
        download_btn = html.Button("Download", id="run-details-download-btn", className="btn-dark mt-3", disabled=True)
    else:
        download_btn = html.Button("Download", id="run-details-download-btn", className="btn-dark mt-3 hidden", disabled=True)
    layout = [
        dbc.Row([
            html.Div(id='hdf5-run-details-url-hash', children="",
                     style={"display": "none"}),

        ]),
        dbc.Row([
            CardCol([
                html.Div([
                    dbc.Card([
                        dbc.CardHeader(dbc.Button(html.H6('Run', id='load-data'), id=f'controle-data', block=True)),
                            dbc.CardBody([
                                html.Div(id='warn-div'),
                                html.Div(id='find-files-warn'),
                                dbc.Row([html.H6('   Runpath'),]),
                                dcc.Input(id='hdf5-run-details-run-path', debounce=True,
                                          value=local_config['default_run_path'], style={"width": "100%"}, placeholder="Enter or paste run path"),
                                html.Div(id='hdf5-run-details-alert'),

                            ])
                    ]),
                    dbc.Card([
                        dbc.CardHeader(dbc.Button(html.H6('Generate h5 Data', id='submit-job-text'),
                                                  id='submit-job-header-button', block=True)),
                        dbc.CardBody([
                            html.Div([
                                html.H6('Login'),
                                dcc.Input(id='hdf5-ssh-username', placeholder='Username', debounce=True,
                                          style={'width': '35%'}),
                                dcc.Input(id='hdf5-ssh-password', placeholder='Password', debounce=True,
                                          style={'width': '35%', 'margin-left': '5px', 'margin-right': '5px'},
                                          type='password'),
                                html.Button('Login', id='hdf5-test-conn',
                                            className='btn-primary', disabled=False, style={'margin-left': '10px'}),
                            ], id='hdf5-login-block'),
                            html.Div([
                                html.H6('HDF5 Workflow Options'),
                                dcc.Input(id='hdf5-lane', placeholder='lane', debounce=True,
                                          style={'width': '35%'}),
                                dcc.Input(id='hdf5-tile', placeholder='Tile', debounce=True,
                                          style={'width': '35%', 'margin-left': '5px', 'margin-right': '5px',
                                                 'margin-top': '5px'}),
                                dcc.Input(id='hdf5-ref', placeholder='genome reference(optional)', debounce=True,
                                          style={'width': '35%', 'margin-top': '5px'}),
                                dcc.Input(id='hdf5-output', placeholder='Output Directory', debounce=True,
                                          style={'width': '35%', 'margin-left': '5px', 'margin-right': '5px',
                                                 'margin-top': '5px'}),
                                dcc.RadioItems(id='hdf5-use-bwa', options=[
                                    {'label': 'Use BWA Aligner', 'value': 'bwa'}
                                ], style={'margin-top': '5px'}, value='bwa')
                            ]),
                            html.H6("Generate HDF5 Files"),
                            html.Br(),
                            dbc.Alert("Please login first.", color="secondary",
                                      id="hdf5-login-prompt"),
                            html.Button("Submit", id="hdf5-submit-btn",
                                        className="btn-primary mr-3", disabled=True),
                            html.Button("Test", id="hdf5-test-btn",
                                        className="btn-warning mr-3", disabled=True),
                            html.Button("My Jobs", id="hdf5-list-jobs-btn",
                                        className="btn-info", disabled=True),

                            dcc.Loading([html.Div(id="hdf5-alerts"),
                                         html.Div(id="hdf5-alerts-2"),
                                         html.Div(id="hdf5-alerts-3"),
                                         html.Div(id='hdf5-alerts-jobs'),
                                         ], type='circle', className="mt-3")
                        ])
                    ]),
                    dbc.Card([
                        dbc.CardHeader(dbc.Button(html.H6('Select/Filter Data', id='filter-data-text'), id=f'filter-data-button', block=True)),
                        dbc.CardBody([
                            dbc.Row([
                                html.Div([
                                    html.H6('Bam h5 File :'),
                                    dcc.Dropdown(
                                        options=[{'label': 'None', 'value': 'None'}],
                                        id='bam-h5-select', disabled=False),
                                ], style={"width": "50%"}),
                                dbc.Button("Load Data", id="load-run-metrics", color='primary'),
                            ]),
                            dbc.Row([
                                html.Div([
                                    html.H6('Tile h5 File:'),
                                    dcc.Input(id='file-select',
                                                 placeholder='choose a file', disabled=True, readOnly=True,),
                                    ]),
                                ]),
                            dbc.Row([
                                dbc.Col([
                                    html.H6('Chastity filter'),
                                    dcc.Input(id='chastity_filter', type='number', debounce=True,
                                              value=None, placeholder=None)
                                    ]),
                                dbc.Col([
                                    html.H6('mapQ filter'),
                                    dcc.Input(id='mapQ_filter', type='number', debounce=True,
                                              value=None, placeholder=None)
                                    ]),
                                dbc.Col([
                                    html.H6('Samples'),
                                    dcc.Input(id='n-samples', type='number', debounce=True,
                                              value=2000, placeholder=None)
                                ]),
                            ]),
                            html.Div([dcc.Slider(
                                        id='cycle-slider',
                                        min=1,
                                        max=518,
                                        step=1,
                                        value=1,
                                        marks={1: '1',
                                               130: '130',
                                               260: '260',
                                               390: '390',
                                               518: '518'}
                                     )], id='slider-div'),
                            html.Div(id='slider-output'),
                            dbc.Row([
                                html.Button('Load data', id='load-button')
                                ])
                            ])
                    ])
                ])
            ], width=4),
            CardCol([
                html.Div(id='out'),
                html.Div([
                    dcc.Graph(id='out-graph')
                ])
            ], width=8)
        ])
    ]
    register_callbacks(app, local_config, datasource, default_metric_values, sample_size)
    return layout

def register_callbacks(app, local_config, datasource, default_metric_values, sample_size):

    @app.callback([Output('hdf5-run-details-run-path', 'value'),
                   Output('hdf5-run-details-alert', 'children'),
                   Output('hdf5-run-details-url-hash', 'children'),
                   ],
                  [Input('url', 'hash')],
                  prevent_initial_call=False)
    def populate_run_from_url_hash(url_hashes):
        if url_hashes is not None and len(url_hashes) > 1:
            data = decode_from_b64(url_hashes[1:])
            app.XXXXXXXXX_logger.debug(
                "URL hash is {}, decoded data is {}".format(url_hashes, data))
            if "Runpath" in data:
                runpath = data["Runpath"]
                if "Charts" in data:
                    alert = dbc.Alert("Shared run loaded, toggle `Live Update` to view".format(
                        basename(runpath)), color="info", dismissable=True)
                else:
                    alert = ""
                metric_sets = data.get("Metrics", default_metric_values)
                return (runpath, alert, json.dumps(data))
            else:
                alert = dbc.Alert(
                    "Cannot find run path in hash or hash decoding failed", color="danger", dismissable=True)
                return None, alert, json.dumps(data)
        else:
            return None, "", default_metric_values, ""

    @app.callback(
        [Output('hdf5-alerts', 'children'),
         Output('hdf5-submit-btn', 'disabled'),
         Output('hdf5-test-btn', 'disabled'),
         Output('hdf5-list-jobs-btn', 'disabled'),
         Output('hdf5-login-prompt', 'is_open'),
         Output('hdf5-login-block', 'style')],
        [Input('hdf5-test-conn', 'n_clicks')],
        [State('hdf5-ssh-username', 'value'),
         State('hdf5-ssh-password', 'value')])
    def test_ssh_connection(clicks, user, password):
        app.XXXXXXXXX_logger.debug(
            "Attempting to connect to cluster with user {}...".format(user))
        server = local_config['ssh_server']
        qsub = QSubManager(server, user=user, password=password)
        app.XXXXXXXXX_logger.debug("Connection {}.".format(
            "alive" if qsub.is_connected() else "failed"))

        if qsub.is_connected():
            alert = dbc.Alert("Connection successful to {}@{}.".format(
                user, server), color="success", dismissable=True)
            return alert, False, False, False, False, {"display": "none"}
        else:
            alert = dbc.Alert("Connection failed to {}@{}.".format(
                user, server), color="danger", dismissable=True)
            return alert, True, True, True, True, {}

    @app.callback([Output('hdf5-alerts-2', 'children'),
                   Output('hdf5-alerts-3', 'children')],
                  [
                      Input('hdf5-submit-btn', 'n_clicks'),
                      Input('hdf5-test-btn', 'n_clicks')
                  ],
                  [State('hdf5-ssh-username', 'value'),
                   State('hdf5-ssh-password', 'value'),
                   State('hdf5-run-details-run-path', 'value'),
                   State('hdf5-lane', 'value'),
                   State('hdf5-tile', 'value'),
                   State('hdf5-ref', 'value'),
                   State('hdf5-output', 'value'),
                   State('hdf5-use-bwa', 'value')])
    def submit_click(nclicks, test_clicks, user, password, runpath, lane, tile, ref, output, aligner):
        if runpath is not None:
            app.XXXXXXXXX_logger.debug("submit runpath: {}".format(runpath))
            runpath = runpath.split(";")[0]
        server = local_config['ssh_server']
        if nclicks is not None or test_clicks is not None:
            is_test = (get_triggered_component()[0] == "hdf5-test-btn")
            runpath = app.path_mapper.map_path(runpath, Platform.Linux)
            msg = "runpath: {}".format(
                runpath)
            app.XXXXXXXXX_logger.debug(
                "Submitting cluster job...{} test? {}".format(msg, is_test))
            submit_success, response = submit_hdf5_generation(
                runpath, user, password, lane, tile, output, ref, aligner, is_test)
            if submit_success:
                alert = dbc.Alert(["Submit successful for bam metrics job to {}@{}.".format(
                    user, server), html.Pre(response)], color="success", dismissable=True)
            else:
                alert = dbc.Alert(["Submit failed for bam metrics job to {}@{}.".format(
                    user, server), html.Pre(response)], color="danger", dismissable=True)

            submit_success_cluster, response_cluster = submit_cluster_h5(user, password, runpath, output, is_test)
            if submit_success_cluster:
                alert2 = dbc.Alert(["Submit successful for cluster h5 job to {}@{}.".format(
                    user, server), html.Pre(response_cluster)], color="success", dismissable=True)
            else:
                alert2 = dbc.Alert(["Submit failed for job cluster h5 to {}@{}.".format(
                    user, server), html.Pre(response_cluster)], color="danger", dismissable=True)

            return [alert, alert2]
        else:
            return ["" ""]

    def submit_hdf5_generation(runpath, user, password, lane, tile, output, ref, aligner, is_test=False, server=local_config['ssh_server']):
        qsub = QSubManager(server, user=user, password=password)

        cmd_ops = ['-r ', '-l ', '-t ', '-o ', '-g ']
        cmd_args = []
        for i, v in enumerate([runpath, lane, tile, output, ref]):
            if v is not None:
                cmd_args.append(cmd_ops[i] + v)
        if output is None:
            output = '.'
        cmd_args = ' '.join(cmd_args)
        cmd = 'qsub -l "excl=true" -N flowCellToh5_Telescope -M {}@illumina.com -m beas -o {}' \
              ' -e {} -b y "/illumina/sync/software/groups/rtaTools/bamMetricsHDF5/fctoh5.sh {}"'.format(user, output+'/qsub.log', output+'/qsub_error.log', cmd_args)
        if aligner == 'bwa':
            cmd = cmd + ' -b'

        if is_test:
            response, exit_code = qsub.exec_command(cmd+' -d')
        else:
            response, exit_code = qsub.exec_command(cmd)

        if exit_code == 0:
            return True, response
        else:
            return False, response

    def submit_cluster_h5(user, password, runpath, output, is_test=False, server=local_config['ssh_server']):
        qsub = QSubManager(server, user=user, password=password)
        cmd_ops = ['-r ', '-o ']
        cmd_args = []
        for i, v in enumerate([runpath, output]):
            if v is not None:
                cmd_args.append(cmd_ops[i]+v)
        cmd_args = ' '.join(cmd_args)

        if output is None:
            output = '.'
            h5_output = ''
        else:
            h5_output = f'-o {output}/neighbour.h5'

        cluster_cmd = f'qsub -l "excl=true" -N cluster_h5_XXXXXXXXX -M {user}@illumina.com -m beas' \
                      f' -o {output}/qsub_findClusterNeighbours.log -e {output}/qsub_findClusterNeighbours_err.log -b y' \
                      f' "/illumina/scratch/SystemsAnalysis/code/scripts/findClusterNeighbours.R {h5_output} {runpath}"'

        if is_test:
            exit_code = 0
            response = 'cluster test successful'
        else:
            response, exit_code = qsub.exec_command(cluster_cmd)

        if exit_code == 0:
            return True, response
        else:
            return False, response

    @app.callback(Output('hdf5-alerts-jobs', 'children'),
        [
            Input('hdf5-list-jobs-btn', 'n_clicks')
        ],
        [
            State('hdf5-ssh-username', 'value'),
            State('hdf5-ssh-password', 'value'),
         ])
    def on_list_jobs_click(nclicks, user, password):
        alert = ""
        if nclicks is not None:
            server = local_config['ssh_server']
            qsub = QSubManager(server, user, password)
            try:
                if qsub.is_connected():
                    result = qsub.list_current_jobs(user)
                    if len(result) > 0:
                        table = get_jobs_list_table(result)
                        alert = dbc.Alert([
                            "Qsub jobs for {}@{}.".format(user, server),
                            html.Pre(table.to_string())], color="info", dismissable=True)
                    else:
                        alert = dbc.Alert("No active/queued jobs found for {}@{}.".format(user, server), color="info", dismissable=True)
                else:
                    alert = dbc.Alert("Cannot query for job info to {}@{}".format(
                        user, server), color="danger", dismissable=True)
            except Exception as e:
                alert = dbc.Alert(["Cannot query for job info to {}@{}".format(user, server),
                                  html.Pre(str(e))], color="danger", dismissable=True)
        return alert


    @app.callback(
        Output('slider-output', 'children'),
        [Input('cycle-slider', 'value')])
    def print_cycle(cycle):
        return 'Cycle: {}'.format(cycle)

    @app.callback([Output('bam-h5-select', 'options')],
                  [Input('hdf5-run-details-run-path', 'value')])
    def bam_h5_selection(runpath):
        files = datasource.find_h5_file(runpath)
        bam_h5s = datasource.list_bam_h5s(files, runpath)
        bam_opt = []
        for file in bam_h5s:
            bam_opt.append({'label': file, 'value': file})

        return [bam_opt]

    @app.callback(
        [Output('file-select', 'value'),
         Output('file-select', 'disabled'),
         Output('find-files-warn', 'children')],
        [Input('load-run-metrics', 'n_clicks')],
        [State('hdf5-run-details-run-path', 'value'),
         State('bam-h5-select', 'value')],
        prevent_initial_call=False)
    def populate_files(clicks, runpath, bam_input):
        if runpath is None or len(runpath) == 0:
            raise PreventUpdate()
        disabled = True
        try:
            app.XXXXXXXXX_logger.debug("Populate files: runpath selected is: {}".format(runpath))
            files = datasource.find_h5_file(runpath)
            if len(files) == 0:
                raise ValueError("Did not find any h5 files")
            else:
                files = datasource.find_h5_file(runpath)
                bam = bam_input
                lane = datasource.get_lane(runpath, bam)
                tile = datasource.get_tile(runpath, lane, bam)
                datasource.load_files(bam_file=bam, tile_file=tile, lane=lane, run=runpath)
                bam_h5s = datasource.select_bam_metric_h5(files, runpath)
                disabled = False
        except Exception as e:
            alert = dbc.Alert(
                'could not find HTF5 files', color="danger", dismissable=True)
            return ['no files found', False,  alert]


        return [tile, disabled, dash.no_update]


    # if files have been found, set the length of the cycle slider
    @app.callback(
        [Output(component_id='cycle-slider', component_property='min'),
         Output(component_id='cycle-slider', component_property='max')],
        [Input('file-select', 'value')])
    def update_slider(file):
        max_cycle = datasource.get_h5_dimension() - 1
        return [1, max_cycle]


    @app.callback(
        Output('out-graph', 'figure'),
        Output('warn-div', 'children'),
        [Input('load-button', 'n_clicks')],
        [State('file-select', 'value'),
         State('mapQ_filter', 'value'),
         State('chastity_filter', 'value'),
         State('cycle-slider', 'value'),
         State('hdf5-run-details-run-path', 'value'),
         State('n-samples', 'value')])
    def plot(clicks, file, mapq, chastity, cycle, run, samples):
        alert = None

        # extract data dependent on selected mode
        data, warn, level = datasource.box_plot_data(run=run, n_cycle=cycle, n_samples=samples)

        # if data extraction fails or if sampling to high notify user
        if warn is not None:
            alert = dbc.Alert(
                warn, color=level, dismissable=True)
            if level == 'danger':
                return [data, alert]

        # filter and plot data
        try:
            if chastity is None and mapq is None:
                fig = px.scatter(data, x='ch1', y='ch2',
                                 color='basecall', facet_col='ref')
            else:
                if chastity is None:
                    chastity = -5
                if mapq is None:
                    mapq = 0
                filtered_data = data[(data.Chastity >= chastity) & (data.mapQ >= mapq)]
                if filtered_data.empty:
                    alert = dbc.Alert(
                        f'filtering resulted in no data, check your filtering options', color='warning', dismissable=True)
                    return [{'data': []}, alert]

                fig = px.scatter(filtered_data, x='ch1', y='ch2',
                                 color='basecall', facet_col='ref')
        except Exception as e:
            alert = dbc.Alert(
                f'failed to plot({e})', color='danger', dismissable=True)
            return [{'data': []}, alert]

        if alert is None:
            return [fig, dash.dash.no_update]
        else:
            return [fig, alert]


def create_tab(app, config, skip_layout=False) -> TelescopeModule:
    APP_NAME = "H5 Bam Metrics"
    APP_DESCRIPTION = ""
    APP_URL = "H5-bam-metrics"
    APP_SHOW_IN_NAVIGATION = True
    APP_SUPPORT_URL_HASH = True
    url = config.get("module_url", APP_URL)
    show_in_navigation = config.get("show_in_navigation", APP_SHOW_IN_NAVIGATION)
    dash_layout = create_layout(app, config, url) if not skip_layout else None
    return TelescopeModule(name=APP_NAME,
                           description=APP_DESCRIPTION,
                           url=url,
                           layout=dash_layout,
                           show_in_navigation=show_in_navigation,
                           support_url_hash=APP_SUPPORT_URL_HASH)

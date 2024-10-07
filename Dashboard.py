#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 13:28:51 2022

@author: grominou
"""

import os
from dash import Dash, dcc, html, Input, Output, State, dash_table, callback_context
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.graph_objects as go
import statsmodels.api as sm
import numpy as np
# import scipy

path = os.getcwd()+"/datasheets"
print(path)


# =============================================================================
# Initialize dash app
# =============================================================================

app = Dash(__name__)
#server = app.server
app.title = "Oxygen isotope equations" #Assigning title to be displayed on tab

# =============================================================================
# Initialize variables
# =============================================================================

color = ['blue', 'red', 'yellow', 'orange', 'green', 'purple', 'grey', 'peru', 'lightblue', 'lightcoral', 'olive']
alpha = .05
predict_placeholder = ["Compute OLS first", "Input phosphate value"]
predicted_value = ""
predicted_point = []

# =============================================================================
# Data loading and cleaning
# =============================================================================

content = os.listdir(path)
xlsx_files = [f for f in content if f[-4:] == "xlsx"]
data_table = pd.DataFrame()

for f in xlsx_files:
    data = pd.read_excel(path+"/"+f)
    data_table = pd.concat([data_table, data], ignore_index = True)
    
data_table = data_table[data_table['Age'].isin(['present-day'])].copy()
data_table.dropna(subset = ["d18Op"], inplace=True)
data_table.drop_duplicates(subset=['Taxon', 'd18Op', 'd18Owmoy'], inplace=True)
data_table['mergedd18Owmoy'] = data_table["d18Owmoy"].fillna(data_table["newd18Owmoy"])

animals = sorted(list(data_table["Animal"].unique()))
ecology = sorted(list(data_table["Environment"].unique()))
water = "mergedd18Owmoy"
marker_color = {}
for i in range(len(animals)):
    marker_color[animals[i]] = color[i]
error = []
ols=[]
pred_env=[]
num = data_table[water].count()
slope_intercept = np.array([0,0])
std_err = np.array([0,0])
equation = 'δ',html.Sup('18'),'O',html.Sub('w'),' = {:.3f}±{:.3f} * δ'.format(slope_intercept[1], std_err[1]),html.Sup('18'),'O',html.Sub('p'),' - {:.3f}±{:.3f}'.format(abs(slope_intercept[0]), std_err[0])


# =============================================================================
# Functions definition for table and graph rendering, ols computation and data prediction
# =============================================================================

def update_dataset(animal, ecology, watersource):
    working_table = data_table[data_table['Animal'].isin(animal)].copy()
    working_table = working_table[working_table['Environment'].isin(ecology)].copy()
    working_table.dropna(subset = [watersource], inplace=True)
    working_table = working_table.sort_values(by=['d18Op']).copy()
    
    return working_table

def create_data_table(working_table, water):
    dat = working_table.to_dict('records')
    
    col = [
        {"name": "Taxon",
         "id": "Taxon",
         "type": "text"
        },
        {"name": "δ¹⁸Oₚ (‰)",
         "id": "d18Op",
         "type": "numeric",
         "format": {"specifier": ",.1f"}
        },
        {"name": "δ¹⁸Oᵥᵥ (‰)",
         "id": water,
         "type": "numeric",
         "format": {"specifier": ",.1f"}
        },
        {"name": "Reference",
         "id": "refShort",
         "type": "text",
        }
        ]
        
    return dat, col
    
def create_scatter(working_table, water, error, ols, pred_env, predicted_point):
        
    fig = go.Figure()
    if error:
        fig.add_trace(go.Scatter(
            x=working_table["d18Op"], 
            y=working_table['newd18Owmoy'], 
            mode='markers', 
            showlegend=False,
            marker=dict(
                color='black',
                size=2
                ),
            error_y=dict(
                type='data',
                symmetric=False,
                array=[
                    [abs(i-j)] for i,j in zip(working_table['newd18Owmoy'], 
                    working_table["newd18Owmax"])
                    ], 
                arrayminus=[
                    [abs(i-j)] for i,j in zip(working_table['newd18Owmoy'], 
                    working_table["newd18Owmin"])
                    ]
                )
            ))
    
    if ols:
        if 'bestfit' in working_table.columns:
            fig.add_trace(go.Scatter(
                name='line of best fit', 
                x=working_table["d18Op"], 
                y=working_table['bestfit'], 
                mode='lines', 
                line_color='red'
                ))
            fig.add_trace(go.Scatter(
                x=working_table["d18Op"], 
                y=working_table['mean_ci_upper'], 
                fill=None, mode='lines', 
                line_color='lightcoral', 
                showlegend=False
                ))
            fig.add_trace(go.Scatter(
                x=working_table["d18Op"], 
                y=working_table['mean_ci_lower'], 
                fill='tonexty', 
                mode='lines', 
                line_color='lightcoral', 
                showlegend=False
                ))       
            if pred_env:
                fig.add_trace(go.Scatter(
                    x=working_table["d18Op"], 
                    y=working_table['obs_ci_upper'], 
                    fill=None, mode='lines', 
                    line_color='lightblue', 
                    showlegend=False
                    ))
                fig.add_trace(go.Scatter(
                    x=working_table["d18Op"], 
                    y=working_table['obs_ci_lower'], 
                    fill='tonexty', 
                    mode='lines', 
                    line_color='lightblue', 
                    showlegend=False
                    ))
            
    for animal, group in working_table.groupby('Animal'):
        fig.add_trace(go.Scatter(
            x=group["d18Op"], 
            y=group[water],
            mode='markers', 
            name=animal,
            marker=dict(
                color=marker_color[animal],
                size=10,
                line=dict(
                    color='black',
                    width=1
                    )
                ),
            customdata=group["Taxon"],
            showlegend=True,
            hovertemplate='<i>%{customdata}</i><br>δ<sup>18</sup>O<sub>p</sub>: %{x:.1f}<br>δ<sup>18</sup>O<sub>w</sub>: %{y:.1f}<extra></extra>'
            ))
    if predicted_point:
        fig.add_trace(go.Scatter(
            x=[predicted_point[0]], 
            y=[predicted_point[1]],
            error_y=dict(
                type='data',
                array=[predicted_point[2]],
                visible=True
                ),
            mode='markers', 
            marker=dict(
                color='black',
                size=10,
                line=dict(
                    color='black',
                    width=1
                    )
                ),
            showlegend=True,
            hovertemplate='<i>%{customdata}</i><br>δ<sup>18</sup>O<sub>p</sub>: %{x:.1f}<br>δ<sup>18</sup>O<sub>w</sub>: %{y:.1f}<extra></extra>'
            ))
    fig.update_xaxes(title_text = "δ<sup>18</sup>O<sub>p</sub> (‰ V-SMOW)")
    fig.update_yaxes(title_text = "δ<sup>18</sup>O<sub>w</sub> (‰ V-SMOW)")
        
    return fig


# =============================================================================
# Initialize tabledata and graph
# =============================================================================

working_table = update_dataset(animals, ecology, water)

dat, col = create_data_table(working_table, water)

fig = create_scatter(working_table, water, error, ols, pred_env, predicted_point)


# =============================================================================
# creating checklist sections
# =============================================================================
 
checklist_boxes = html.Div(
    id = 'checklist_boxes',
    children = [
        html.Div(
            id = "bloc_animals",
            children=[
                html.H2('Select animals'),
                    
                html.Button(
                    'select all / none',
                    id = "all_or_none1",
                    n_clicks=0
                    ),
                    
                dcc.Checklist(
                    id = "selected_animals",
                    options=[{"label": " " + i, "value": i} for i in animals],
                    value = animals,
                    labelStyle={'display': 'block'},
                    )
                ],
            ),
        
        html.Div(
            id = "bloc_ecology",
            children=[
                html.H2('Select ecology'),
                    
                html.Button(
                    'select all / none',
                    id = "all_or_none2",
                    n_clicks=0
                    ),
                    
                dcc.Checklist(
                    id = "selected_ecology",
                    options=[{"label": " " + i, "value": i} for i in ecology],
                    value = ecology,
                    labelStyle={'display': 'block'},
                    )
                ],
            ),
        
        html.Div(
            id = "bloc_water",
            children=[
                html.H2('Select water data source'),
                dcc.RadioItems(
                    id='selected_water',
                    options=[
                        {"label": " Published + New", "value": "pubnew"},
                        {"label": " Published data", "value": "original"},
                        {"label": " New data", "value": "new"},
                        ],
                    value="pubnew",
                    labelStyle={'display': 'block'}
                    )
                ],
            ),  
    ])


# =============================================================================
# creating data table section
# =============================================================================

datatable = html.Div(
    id = "bloc_datatable",
    children=[
        html.H2('data table'),
        dash_table.DataTable(
            id='table',
            columns=col,
            data=dat,
            row_selectable="multi",
            selected_rows=list(range(0, len(working_table))),
            style_cell_conditional=[
                    {
                    'if': {'column_id': "Taxon"},
                    'textAlign': 'left',
                    'font-style': 'italic'
                    }
                ],
            style_cell={'font-family': 'sans-serif'},
            style_header={'fontWeight': 'bold', 'font-style': 'normal'},
            style_data={'minWidth': '20px', 'width': 'auto', 'maxWidth': '500px'},
            page_size=1000,
            style_table={'height': '55vh', 'width': '30vw', 'overflowY': 'auto'}
            ),
        html.P(
            'N = ' + str(num),
            id = "count",
            ),
        ])

# =============================================================================
# creating graph section
# =============================================================================

graph = html.Div(children = [
    html.Div(
        id = "bloc_plot",
        children=[
            html.H2('Phosphate-water plot'),
            
            dcc.Graph(
                id='graph',
                figure=fig,
                style = {'width': '47vw', 'height': '60vh'
                    },
                ),
            dcc.Checklist(
                id="error_bars",
                options=[{"label": " Environmental min/max", "value": "minmax"}],
                value = [],
                )
            ]),
        ])

# =============================================================================
# creating line of best fit and prediction sections
# =============================================================================

bestfit = html.Div(children = [
    html.Div(
        id = "bloc_stat",
        children=[
            html.H2('Line of best fit'),
            
            dcc.Checklist(
                id = "ols",
                options=[{"label": ' Compute OLS', "value": "ols"}],
                value=[]
                ),
            
            dcc.Checklist(
                id = "pred_env",
                options=[{"label": ' Predict. env.', "value": "pred_env"}],
                value=[]
                ),
            
            html.P(
                equation,
                id = "ols_equation",
                ),
            ]),
        ])

prediction = html.Div(children = [
    html.Div(
        id = "bloc_predict",
        children=[
            html.H2('Water value Prediction'),

            html.P(
                ('δ',html.Sup('18'),'O',html.Sub('p'),'value (‰ V-SMOW): ')
                ),
            
            dcc.Input(
                id = "prediction",
                type = "number",
                disabled = True,
                placeholder = predict_placeholder[0],
                debounce=True,
                ),

            html.P(
                ('Predicted δ',html.Sup('18'),'O',html.Sub('w'),'value (‰ V-SMOW): ')
                ),
            
            html.P(
                predicted_value,
                id = "predicted_val",
                ),
            ]),
        ])
 

# =============================================================================
# adding the layout to the dash app
# =============================================================================

app.layout = html.Div(
    id = "main_div", 
    children =[checklist_boxes, datatable, graph, bestfit, prediction]
    )


# =============================================================================
# Application callbacks
# =============================================================================

@app.callback(
    Output('table', 'data'),
    Output('table', 'columns'),
    Output('graph', 'figure'),
    Output('count', 'children'),
    Output('ols_equation', 'children'),
    Output('prediction', 'placeholder'),
    Output('prediction', 'disabled'),
    Output('predicted_val', 'children'),
    Input('selected_animals', 'value'),
    Input('selected_ecology', 'value'),
    Input('selected_water', 'value'),
    Input('error_bars', 'value'),
    Input("ols", "value"),
    Input("pred_env", "value"),
    Input("prediction", "value")
    )

def update_tableandfig(animal_value, ecology_value, water_value, err, ols, pred_env, prediction):
        
    if water_value == 'new':
        water = 'newd18Owmoy'
    elif water_value == 'original':
        water = 'd18Owmoy'         
    else:
        water = 'mergedd18Owmoy'
    
    working_table = update_dataset(animal_value, ecology_value, water)
    
    if ols and len(working_table) > 1:
        model = sm.OLS(working_table[water],sm.add_constant(working_table['d18Op']))
        model_results = model.fit()
        predictions = model_results.get_prediction().summary_frame(alpha)
        working_table['bestfit'] = model_results.fittedvalues
        working_table['mean_ci_lower'] = predictions['mean_ci_lower']
        working_table['mean_ci_upper'] = predictions['mean_ci_upper']
        slope_intercept = model_results.params
        std_err = model_results.bse
        equation = 'δ',html.Sup('18'),'O',html.Sub('w'),' = {:.3f}±{:.3f} * δ'.format(slope_intercept.iloc[1], std_err.iloc[1]),html.Sup('18'),'O',html.Sub('p'),' - {:.3f}±{:.3f}'.format(abs(slope_intercept.iloc[0]), std_err.iloc[0])
        placeholder = predict_placeholder[1]
        disabled = False
        
        if pred_env:
            working_table['obs_ci_lower'] = predictions['obs_ci_lower']
            working_table['obs_ci_upper'] = predictions['obs_ci_upper']
    
    else:
        equation = ''
        placeholder = predict_placeholder[0]
        disabled = True

    dataset, column = create_data_table(working_table, water)

    if isinstance(prediction, float) or isinstance(prediction, int):
        dataForPrediction = np.array([[1, prediction]])
        predict = model_results.get_prediction(dataForPrediction)
        pred_summary = predict.summary_frame(alpha)
        ypred = pred_summary['mean'][0]
        yprederror = pred_summary['mean'][0]-pred_summary['obs_ci_lower'][0]
        predicted_value = '{:.1f} ± {:.1f}'.format(ypred, yprederror)
        predicted_point = [prediction, ypred, yprederror]

    else:
        predicted_value = ""
        predicted_point = []


    fig = create_scatter(working_table, water, err, ols, pred_env, predicted_point)
    num = working_table[water].count()
    count = 'N = ' + str(num)
    
    return dataset, column, fig, count, equation, placeholder, disabled, predicted_value

@app.callback(
    Output("selected_animals", "value"),
    Input("all_or_none1", "n_clicks"),
    State("selected_animals", "options")
    )
   
def all_or_none1(n_clicks, options):
    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate
    else:
        trigged_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
        if trigged_id == 'all_or_none1':
            if n_clicks == 0: ## Do not update in the beginning
                raise PreventUpdate
    
            if n_clicks % 2 == 0: ## Select all options on even clicks
                return [i['value'] for i in options]
            else: ## Clear all options on odd clicks
                return []
        else:
            raise PreventUpdate

@app.callback(
    Output("selected_ecology", "value"),
    Input("all_or_none2", "n_clicks"),
    State("selected_ecology", "options")
    )
   
def all_or_none2(n_clicks, options):
    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate
    else:
        trigged_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
        if trigged_id == 'all_or_none2':
            if n_clicks == 0: ## Do not update in the beginning
                raise PreventUpdate
    
            if n_clicks % 2 == 0: ## Clear all options on even clicks
                return [i['value'] for i in options]
            else: ## Select all options on odd clicks
                return []
        else:
            raise PreventUpdate


# =============================================================================
# Run application
# =============================================================================

if __name__ == '__main__':
    app.run_server(debug=False, port = 8080)

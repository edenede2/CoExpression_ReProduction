"""
Visualization Module
Creates plots and visualizations for network analysis and MDC results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from typing import Dict, List, Optional, Tuple
import warnings

class NetworkVisualizer:
    """Class for creating network analysis visualizations."""
    
    def __init__(self):
        """Initialize the visualizer."""
        self.color_palette = px.colors.qualitative.Set3
        
    def plot_power_analysis(self, power_analysis_df: pd.DataFrame) -> go.Figure:
        """
        Plot soft threshold power analysis results.
        
        Args:
            power_analysis_df: Power analysis results DataFrame
            
        Returns:
            Plotly figure
        """
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Scale-free topology fit", "Mean connectivity"),
            x_title="Soft threshold power"
        )
        
        # Scale-free topology R²
        fig.add_trace(
            go.Scatter(
                x=power_analysis_df['power'],
                y=power_analysis_df['rsq'],
                mode='lines+markers',
                name='R²',
                line=dict(color='blue')
            ),
            row=1, col=1
        )
        
        # Add horizontal line at R² = 0.8
        fig.add_hline(y=0.8, line_dash="dash", line_color="red", row=1, col=1)
        
        # Mean connectivity
        fig.add_trace(
            go.Scatter(
                x=power_analysis_df['power'],
                y=power_analysis_df['mean_connectivity'],
                mode='lines+markers',
                name='Mean connectivity',
                line=dict(color='green')
            ),
            row=1, col=2
        )
        
        fig.update_layout(
            title="Soft Threshold Power Selection",
            height=400,
            showlegend=False
        )
        
        fig.update_yaxes(title_text="R²", row=1, col=1)
        fig.update_yaxes(title_text="Mean connectivity", row=1, col=2)
        
        return fig
    
    def plot_sample_pca(self, pca_results: pd.DataFrame, color_by: Optional[str] = None) -> go.Figure:
        """
        Plot PCA of samples with outlier detection.
        
        Args:
            pca_results: PCA results DataFrame
            color_by: Column to color points by
            
        Returns:
            Plotly figure
        """
        if color_by and color_by in pca_results.columns:
            color = pca_results[color_by]
        else:
            color = pca_results['is_outlier'].map({True: 'Outlier', False: 'Normal'})
        
        fig = px.scatter(
            pca_results,
            x='PC1',
            y='PC2',
            color=color,
            title="Sample PCA with Outlier Detection",
            hover_data=['mahalanobis_distance']
        )
        
        fig.update_layout(
            xaxis_title=f"PC1 ({pca_results.columns[0]})",
            yaxis_title=f"PC2 ({pca_results.columns[1]})",
            height=500
        )
        
        return fig
    
    def plot_module_sizes(self, modules: Dict[str, int]) -> go.Figure:
        """
        Plot module size distribution.
        
        Args:
            modules: Gene-to-module mapping
            
        Returns:
            Plotly figure
        """
        module_counts = pd.Series(modules).value_counts().sort_index()
        
        # Remove grey module (0) if present
        if 0 in module_counts.index:
            grey_count = module_counts[0]
            module_counts = module_counts.drop(0)
        else:
            grey_count = 0
        
        fig = go.Figure(data=[
            go.Bar(
                x=module_counts.index,
                y=module_counts.values,
                name='Module sizes'
            )
        ])
        
        fig.update_layout(
            title=f"Module Size Distribution (Grey module: {grey_count} genes)",
            xaxis_title="Module ID",
            yaxis_title="Number of genes",
            height=400
        )
        
        return fig
    
    def plot_correlation_heatmap(self, 
                                correlation_matrix: pd.DataFrame,
                                max_genes: int = 100) -> go.Figure:
        """
        Plot gene correlation heatmap.
        
        Args:
            correlation_matrix: Gene correlation matrix
            max_genes: Maximum number of genes to plot
            
        Returns:
            Plotly figure
        """
        # Subsample if too many genes
        if correlation_matrix.shape[0] > max_genes:
            sample_indices = np.random.choice(
                correlation_matrix.index, 
                size=max_genes, 
                replace=False
            )
            plot_matrix = correlation_matrix.loc[sample_indices, sample_indices]
        else:
            plot_matrix = correlation_matrix
        
        fig = go.Figure(data=go.Heatmap(
            z=plot_matrix.values,
            x=plot_matrix.columns,
            y=plot_matrix.index,
            colorscale='RdBu_r',
            zmid=0,
            colorbar=dict(title="Correlation")
        ))
        
        fig.update_layout(
            title=f"Gene Correlation Heatmap ({len(plot_matrix)} genes)",
            height=600,
            xaxis_showticklabels=False,
            yaxis_showticklabels=False
        )
        
        return fig
    
    def plot_module_eigengenes(self, 
                               eigengenes_df: pd.DataFrame,
                               sample_groups: Optional[Dict[str, str]] = None) -> go.Figure:
        """
        Plot module eigengenes across samples.
        
        Args:
            eigengenes_df: Module eigengenes DataFrame
            sample_groups: Dictionary mapping sample IDs to group labels
            
        Returns:
            Plotly figure
        """
        # Prepare data for plotting
        plot_data = eigengenes_df.reset_index().melt(
            id_vars='index',
            var_name='Module',
            value_name='Eigengene'
        )
        plot_data.rename(columns={'index': 'Sample'}, inplace=True)
        
        # Add group information if provided
        if sample_groups:
            plot_data['Group'] = plot_data['Sample'].map(sample_groups)
            color = 'Group'
        else:
            color = None
        
        fig = px.box(
            plot_data,
            x='Module',
            y='Eigengene',
            color=color,
            title="Module Eigengenes Across Samples"
        )
        
        fig.update_layout(
            xaxis_title="Module",
            yaxis_title="Eigengene Value",
            height=500
        )
        
        return fig
    
    def plot_mdc_results(self, mdc_results: pd.DataFrame, 
                         condition_1_name: str = "Condition1",
                         condition_2_name: str = "Condition2") -> go.Figure:
        """
        Plot MDC analysis results.
        
        Args:
            mdc_results: MDC results DataFrame
            condition_1_name: Name of condition 1
            condition_2_name: Name of condition 2
            
        Returns:
            Plotly figure
        """
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "MDC Ratio vs Module Size",
                "Connectivity Comparison",
                "Log MDC Ratio Distribution",
                "Density Difference"
            )
        )
        
        # MDC ratio vs module size
        fig.add_trace(
            go.Scatter(
                x=mdc_results['module_size'],
                y=mdc_results['mdc_ratio'],
                mode='markers',
                name='MDC Ratio',
                text=mdc_results['module_id'],
                hovertemplate='Module %{text}<br>Size: %{x}<br>MDC Ratio: %{y:.3f}'
            ),
            row=1, col=1
        )
        
        # Connectivity comparison
        fig.add_trace(
            go.Scatter(
                x=mdc_results[f'{condition_1_name}_connectivity'],
                y=mdc_results[f'{condition_2_name}_connectivity'],
                mode='markers',
                name='Connectivity',
                text=mdc_results['module_id'],
                hovertemplate='Module %{text}<br>%s: %{x:.3f}<br>%s: %{y:.3f}' % (condition_1_name, condition_2_name)
            ),
            row=1, col=2
        )
        
        # Add diagonal line for connectivity comparison
        max_conn = max(
            mdc_results[f'{condition_1_name}_connectivity'].max(),
            mdc_results[f'{condition_2_name}_connectivity'].max()
        )
        fig.add_trace(
            go.Scatter(
                x=[0, max_conn],
                y=[0, max_conn],
                mode='lines',
                line=dict(dash='dash', color='grey'),
                showlegend=False
            ),
            row=1, col=2
        )
        
        # Log MDC ratio distribution
        fig.add_trace(
            go.Histogram(
                x=mdc_results['log_mdc_ratio'].dropna(),
                name='Log MDC Ratio',
                nbinsx=20
            ),
            row=2, col=1
        )
        
        # Density difference
        fig.add_trace(
            go.Scatter(
                x=mdc_results['module_size'],
                y=mdc_results['density_difference'],
                mode='markers',
                name='Density Diff',
                text=mdc_results['module_id'],
                hovertemplate='Module %{text}<br>Size: %{x}<br>Density Diff: %{y:.3f}'
            ),
            row=2, col=2
        )
        
        # Add horizontal line at y=0 for density difference
        fig.add_hline(y=0, line_dash="dash", line_color="grey", row=2, col=2)
        
        fig.update_layout(
            title="Modular Differential Connectivity Analysis",
            height=800,
            showlegend=False
        )
        
        # Update axis labels
        fig.update_xaxes(title_text="Module Size", row=1, col=1)
        fig.update_yaxes(title_text="MDC Ratio", row=1, col=1)
        
        fig.update_xaxes(title_text=f"{condition_1_name} Connectivity", row=1, col=2)
        fig.update_yaxes(title_text=f"{condition_2_name} Connectivity", row=1, col=2)
        
        fig.update_xaxes(title_text="Log MDC Ratio", row=2, col=1)
        fig.update_yaxes(title_text="Frequency", row=2, col=1)
        
        fig.update_xaxes(title_text="Module Size", row=2, col=2)
        fig.update_yaxes(title_text="Density Difference", row=2, col=2)
        
        return fig
    
    def plot_significant_modules(self, 
                                significant_modules: pd.DataFrame,
                                top_n: int = 10) -> go.Figure:
        """
        Plot top significant modules from MDC analysis.
        
        Args:
            significant_modules: Significant modules DataFrame
            top_n: Number of top modules to show
            
        Returns:
            Plotly figure
        """
        if len(significant_modules) == 0:
            # Create empty plot
            fig = go.Figure()
            fig.add_annotation(
                text="No significant modules found",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            fig.update_layout(title="Significant Modules")
            return fig
        
        # Take top N modules by absolute log ratio
        top_modules = significant_modules.head(top_n).copy()
        
        # Create color based on direction of change
        colors = ['red' if x > 0 else 'blue' for x in top_modules['log_mdc_ratio']]
        
        fig = go.Figure(data=[
            go.Bar(
                x=top_modules['module_id'],
                y=top_modules['log_mdc_ratio'],
                marker_color=colors,
                text=top_modules['p_value_corrected'].round(4),
                textposition='outside',
                hovertemplate='Module %{x}<br>Log MDC Ratio: %{y:.3f}<br>Adj. p-value: %{text}'
            )
        ])
        
        fig.update_layout(
            title=f"Top {len(top_modules)} Significant Modules",
            xaxis_title="Module ID",
            yaxis_title="Log₂ MDC Ratio",
            height=500
        )
        
        # Add horizontal line at y=0
        fig.add_hline(y=0, line_dash="dash", line_color="grey")
        
        return fig
    
    def create_network_plot(self, 
                            adjacency_matrix: pd.DataFrame,
                            modules: Dict[str, int],
                            layout: str = 'spring',
                            max_nodes: int = 100) -> go.Figure:
        """
        Create network visualization.
        
        Args:
            adjacency_matrix: Gene adjacency matrix
            modules: Gene-to-module mapping
            layout: Network layout algorithm
            max_nodes: Maximum number of nodes to plot
            
        Returns:
            Plotly figure
        """
        # Sample nodes if too many
        if adjacency_matrix.shape[0] > max_nodes:
            sample_nodes = np.random.choice(
                adjacency_matrix.index,
                size=max_nodes,
                replace=False
            )
            adj_sample = adjacency_matrix.loc[sample_nodes, sample_nodes]
        else:
            adj_sample = adjacency_matrix
            sample_nodes = adjacency_matrix.index
        
        # Create NetworkX graph
        G = nx.Graph()
        
        # Add nodes
        for node in adj_sample.index:
            module_id = modules.get(node, 0)
            G.add_node(node, module=module_id)
        
        # Add edges (only strong connections)
        threshold = adj_sample.values.flatten()
        threshold = np.percentile(threshold[threshold > 0], 95)  # Top 5% of connections
        
        for i, node1 in enumerate(adj_sample.index):
            for j, node2 in enumerate(adj_sample.columns):
                if i < j and adj_sample.iloc[i, j] > threshold:
                    G.add_edge(node1, node2, weight=adj_sample.iloc[i, j])
        
        # Calculate layout
        if layout == 'spring':
            pos = nx.spring_layout(G, k=1, iterations=50)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        else:
            pos = nx.random_layout(G)
        
        # Prepare node traces by module
        node_traces = []
        unique_modules = set(modules.get(node, 0) for node in sample_nodes)
        
        for module_id in unique_modules:
            module_nodes = [node for node in sample_nodes if modules.get(node, 0) == module_id]
            
            if len(module_nodes) == 0:
                continue
            
            node_x = [pos[node][0] for node in module_nodes if node in pos]
            node_y = [pos[node][1] for node in module_nodes if node in pos]
            
            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers',
                name=f'Module {module_id}',
                marker=dict(
                    size=8,
                    color=self.color_palette[module_id % len(self.color_palette)]
                ),
                text=module_nodes,
                hovertemplate='%{text}<br>Module: ' + str(module_id)
            )
            node_traces.append(node_trace)
        
        # Prepare edge trace
        edge_x = []
        edge_y = []
        
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='lightgrey'),
            hoverinfo='none',
            mode='lines',
            showlegend=False
        )
        
        # Create figure
        fig = go.Figure(data=[edge_trace] + node_traces)
        
        fig.update_layout(
            title=f"Gene Coexpression Network ({len(sample_nodes)} genes)",
            showlegend=True,
            hovermode='closest',
            margin=dict(b=20,l=5,r=5,t=40),
            annotations=[
                dict(
                    text=f"Showing top {len(G.edges())} connections",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002,
                    xanchor='left', yanchor='bottom',
                    font=dict(size=12)
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=600
        )
        
        return fig

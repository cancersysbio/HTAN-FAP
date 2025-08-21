import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize
import seaborn as sns

# Column name mappings - set these to match your input data format
COLUMN_NAMES = {
    'alt_depth': 't_alt_count',  # Standard name = actual name in data
    'total_depth': 't_depth'
}

# Set up plotting
plt.style.use('seaborn-v0_8')
sns.set_context('notebook')

def truncated_binom_pmf(k, n, p, truncation_threshold=None):
    """Calculate the probability mass function of a truncated binomial distribution.
    
    Parameters:
    -----------
    k : array-like
        Number of successes (alt allele counts)
    n : array-like
        Number of trials (total depth)
    p : float
        Probability of success (VAF)
    truncation_threshold : float, optional
        Minimum VAF threshold for truncation. If None, uses minimum observed VAF.
    
    Returns:
    --------
    float
        Probability mass
    """
    # If no truncation threshold provided, use minimum observed VAF
    if truncation_threshold is None:
        truncation_threshold = np.min(k / n)
    
    # Calculate observed VAFs
    vafs = k / n
    
    # Zero out probabilities for values below threshold
    mask = vafs >= truncation_threshold
    
    # Calculate raw binomial probabilities
    probs = stats.binom.pmf(k, n, p)
    
    # Zero out probabilities below threshold
    probs[~mask] = 0
    
    # Calculate truncation probability (probability of being below threshold)
    trunc_prob = stats.binom.cdf(np.ceil(n * truncation_threshold) - 1, n, p)
    
    # Return normalized probabilities
    return probs / (1 - trunc_prob)

def mixture_log_likelihood(params, k, n, n_components, truncation_threshold=None):
    """Calculate the log likelihood of the mixture model.
    
    Parameters:
    -----------
    params : array-like
        Model parameters (weights and VAFs)
    k : array-like
        Alt allele counts
    n : array-like
        Total depths
    n_components : int
        Number of mixture components
    truncation_threshold : float, optional
        Minimum VAF threshold for truncation. If None, uses minimum observed VAF.
    
    Returns:
    --------
    float
        Negative log likelihood (for minimization)
    """
    # Split parameters into weights and VAFs
    weights = params[:n_components]
    vafs = params[n_components:]
    
    # Normalize weights using the log-sum-exp trick to avoid overflow
    weights_exp = weights - np.max(weights)  # Subtract max for numerical stability
    weights = np.exp(weights_exp) / np.sum(np.exp(weights_exp))
    
    # Calculate mixture probabilities
    mixture_probs = np.zeros(len(k))
    for i in range(n_components):
        comp_probs = truncated_binom_pmf(k, n, vafs[i], truncation_threshold=truncation_threshold)
        # Add to mixture with weight
        mixture_probs += weights[i] * comp_probs
    
    # Return negative log likelihood with a small epsilon to avoid log(0)
    return -np.sum(np.log(mixture_probs + 1e-10))

def fit_mixture_model(k, n, n_components, n_restarts=10, truncation_threshold=None):
    """Fit the mixture model with specified number of components.
    
    Parameters:
    -----------
    k : array-like
        Alt allele counts
    n : array-like
        Total depths
    n_components : int
        Number of mixture components
    n_restarts : int
        Number of random restarts for optimization
    truncation_threshold : float, optional
        Minimum VAF threshold for truncation. If None, uses minimum observed VAF.
    
    Returns:
    --------
    dict
        Fitted model parameters and statistics
    """
    best_ll = float('inf')
    best_params = None
    
    # Calculate minimum observed VAF for initialization
    min_vaf = np.min(k / n)
    
    for _ in range(n_restarts):
        # Initialize parameters randomly
        init_weights = np.random.dirichlet(np.ones(n_components))
        # Initialize VAFs between min_vaf and 0.5
        init_vafs = np.random.uniform(min_vaf, 0.5, n_components)
        init_params = np.concatenate([np.log(init_weights), init_vafs])
        
        # Optimize parameters
        bounds = [(None, None)] * n_components + [(min_vaf, 0.5)] * n_components
        result = minimize(mixture_log_likelihood, init_params, 
                         args=(k, n, n_components, truncation_threshold),
                         bounds=bounds,
                         method='L-BFGS-B')
        
        if result.fun < best_ll:
            best_ll = result.fun
            best_params = result.x
    
    # Extract final parameters
    weights = np.exp(best_params[:n_components])
    weights = weights / np.sum(weights)
    vafs = best_params[n_components:]
    
    # Calculate AIC and BIC
    n_params = 2 * n_components - 1  # weights sum to 1, so one less parameter
    aic = 2 * best_ll + 2 * n_params
    bic = 2 * best_ll + n_params * np.log(len(k))
    
    return {
        'weights': weights,
        'vafs': vafs,
        'log_likelihood': -best_ll,
        'aic': aic,
        'bic': bic,
        'n_components': n_components,
        'min_vaf': min_vaf if truncation_threshold is None else truncation_threshold
    }

def plot_mixture_model(k, n, model, ax=None):
    """Plot the fitted mixture model and data distribution."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Plot histogram of VAFs
    vafs = k / n
    hist, bin_edges = np.histogram(vafs, bins=50, density=True)
    ax.hist(vafs, bins=50, density=True, alpha=0.3, label='Observed VAFs')
    
    # Get the max height of the histogram for scaling
    hist_max = np.max(hist)
    
    # Plot mixture components
    x = np.linspace(0, 1, 1000)
    n_mean = n.mean()
    
    # Calculate the total PMF to find its maximum value
    y_total = np.zeros_like(x)
    for i in range(model['n_components']):
        y_comp = model['weights'][i] * truncated_binom_pmf(np.round(x * n_mean), np.round(n_mean), model['vafs'][i], truncation_threshold=model['min_vaf'])
        y_total += y_comp
    
    # Find the max of the mixture model for scaling
    pmf_max = np.max(y_total)
    
    # Scale factor to match heights
    scale_factor = hist_max / pmf_max if pmf_max > 0 else 1
    
    # Calculate variants per component if available
    variants_per_component = None
    if 'variants_per_component' in model:
        variants_per_component = model['variants_per_component']
    
    # Plot scaled mixture components
    for i in range(model['n_components']):
        y = model['weights'][i] * truncated_binom_pmf(np.round(x * n_mean), np.round(n_mean), model['vafs'][i], truncation_threshold=model['min_vaf'])
        y_scaled = y * scale_factor
        
        # Create label with variant count if available
        label = f'Component {i+1} (VAF={model["vafs"][i]:.3f}, weight={model["weights"][i]:.3f}'
        if variants_per_component is not None:
            n_variants = int(round(np.nan_to_num(variants_per_component[i], nan=0)))
            label += f', n={n_variants}'
        label += ')'
        
        ax.plot(x, y_scaled, label=label)
    
    # Plot scaled total mixture
    ax.plot(x, y_total * scale_factor, 'k--', label='Total mixture')
    
    ax.set_xlabel('VAF')
    ax.set_ylabel('Density')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15))
    plt.subplots_adjust(bottom=0.2)
    return ax

def is_clonal(model, threshold=0.3):
    """Determine if a sample is clonal based on the VAF of the dominant component.
    
    A sample is considered clonal if the VAF of the dominant component (component with highest weight)
    exceeds the specified threshold.
    
    Parameters:
    -----------
    model : dict
        Fitted mixture model parameters from fit_mixture_model
    threshold : float, optional
        VAF threshold for clonality (default: 0.3)
        
    Returns:
    --------
    bool
        True if the sample is clonal, False otherwise
    float
        VAF of the dominant component
    """
    # Find the index of the component with highest weight
    dominant_idx = np.argmax(model['weights'])
    dominant_vaf = model['vafs'][dominant_idx]
    return dominant_vaf >= threshold, dominant_vaf

def estimate_variants_per_component(k, n, model):
    """Estimate the number of variants belonging to each component of the mixture model.
    
    Parameters:
    -----------
    k : array-like
        Alt allele counts
    n : array-like
        Total depths
    model : dict
        Fitted mixture model parameters from fit_mixture_model
        
    Returns:
    --------
    array-like
        Number of variants estimated to belong to each component
    """
    # Calculate probabilities for each variant belonging to each component
    probs = np.zeros((len(k), model['n_components']))
    for i in range(model['n_components']):
        probs[:, i] = model['weights'][i] * truncated_binom_pmf(k, n, model['vafs'][i], truncation_threshold=model['min_vaf'])
    
    # Normalize probabilities for each variant
    probs = probs / probs.sum(axis=1, keepdims=True)
    
    # Sum probabilities to get expected number of variants per component
    variants_per_component = probs.sum(axis=0)
    
    return variants_per_component

def fit_models_for_sample(data, clonality_threshold=0.3):
    """Fit mixture models for a single sample and return model parameters.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame containing alt_depth and total_depth columns for a single sample
    clonality_threshold : float, optional
        VAF threshold for determining clonality (default: 0.3)
        
    Returns:
    --------
    dict
        Dictionary containing model parameters and statistics
    """
    k = data[COLUMN_NAMES['alt_depth']].to_numpy()
    n = data[COLUMN_NAMES['total_depth']].to_numpy()
    
    # Fit models with different numbers of components
    truncated_models = {}
    non_truncated_models = {}
    for n_components in [1, 2, 3]:
        truncated_models[n_components] = fit_mixture_model(k, n, n_components)
        non_truncated_models[n_components] = fit_mixture_model(k, n, n_components, truncation_threshold=0)
    
    # Find best models based on AIC
    truncated_aics = {n: m['aic'] for n, m in truncated_models.items()}
    non_truncated_aics = {n: m['aic'] for n, m in non_truncated_models.items()}
    
    best_truncated_n = min(truncated_aics, key=truncated_aics.get)
    best_non_truncated_n = min(non_truncated_aics, key=non_truncated_aics.get)
    
    best_truncated = truncated_models[best_truncated_n]
    best_non_truncated = non_truncated_models[best_non_truncated_n]
    
    # Get clonality information
    is_clonal_truncated, dominant_vaf_truncated = is_clonal(best_truncated, clonality_threshold)
    is_clonal_non_truncated, dominant_vaf_non_truncated = is_clonal(best_non_truncated, clonality_threshold)
    
    # Estimate variants per component
    truncated_variants = estimate_variants_per_component(k, n, best_truncated)
    non_truncated_variants = estimate_variants_per_component(k, n, best_non_truncated)
    
    # Calculate mean depth
    mean_depth = np.mean(n)
    
    # Store results
    return {
        'truncated_model': best_truncated,
        'non_truncated_model': best_non_truncated,
        'truncated_n_components': best_truncated_n,
        'non_truncated_n_components': best_non_truncated_n,
        'truncated_aic': best_truncated['aic'],
        'non_truncated_aic': best_non_truncated['aic'],
        'is_clonal_truncated': is_clonal_truncated,
        'is_clonal_non_truncated': is_clonal_non_truncated,
        'dominant_vaf_truncated': dominant_vaf_truncated,
        'dominant_vaf_non_truncated': dominant_vaf_non_truncated,
        'truncated_variants_per_component': truncated_variants,
        'non_truncated_variants_per_component': non_truncated_variants,
        'total_variants': len(k),
        'mean_depth': mean_depth
    }

def plot_clones(data, model_params, name, clonality_threshold=0.3):
    """Plot histogram of VAFs with fitted mixture model components using pre-fitted models.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame containing alt_depth and total_depth columns
    model_params : dict
        Dictionary containing model parameters and statistics from fit_models_for_sample
    name : str
        Title for the plot
    clonality_threshold : float, optional
        VAF threshold for determining clonality (default: 0.3)
    """
    k = data[COLUMN_NAMES['alt_depth']].to_numpy()
    n = data[COLUMN_NAMES['total_depth']].to_numpy()
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Add variant counts to models
    truncated_model = model_params['truncated_model'].copy()
    non_truncated_model = model_params['non_truncated_model'].copy()
    truncated_model['variants_per_component'] = model_params['truncated_variants_per_component']
    non_truncated_model['variants_per_component'] = model_params['non_truncated_variants_per_component']
    
    # Plot truncated model
    plot_mixture_model(k, n, truncated_model, ax=ax1)
    clonality_status_truncated = "Clonal" if model_params['is_clonal_truncated'] else "Not clonal"
    ax1.set_title(f"Truncated Model\n{clonality_status_truncated}, dominant VAF: {model_params['dominant_vaf_truncated']:.3f}")
    
    # Plot non-truncated model
    plot_mixture_model(k, n, non_truncated_model, ax=ax2)
    clonality_status_non_truncated = "Clonal" if model_params['is_clonal_non_truncated'] else "Not clonal"
    ax2.set_title(f"Non-truncated Model\n{clonality_status_non_truncated}, dominant VAF: {model_params['dominant_vaf_non_truncated']:.3f}")
    
    # Add sample name, total variants, and mean depth to the figure
    fig.suptitle(f"{name} (Total variants: {model_params['total_variants']}, Mean depth: {model_params['mean_depth']:.1f}x)", y=1.05)
    
    plt.tight_layout()
    plt.show() 
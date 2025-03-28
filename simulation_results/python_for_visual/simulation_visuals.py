import pandas as pd
import matplotlib.pyplot as plt

# Construct the data frame from the simulation results
data = [
    # n = 100, sim_data = "po"
    {"n": 100, "sim_data": "po", "fit_model": "po", "MSE": 0.2178865, "MAD": 0.3829917, "TV": 0.4650997},
    {"n": 100, "sim_data": "po", "fit_model": "ph", "MSE": 0.2246444, "MAD": 0.3880761, "TV": 0.4775318},
    {"n": 100, "sim_data": "po", "fit_model": "np", "MSE": 0.2048508, "MAD": 0.3622654, "TV": 0.4422477},
    # n = 100, sim_data = "ph"
    {"n": 100, "sim_data": "ph", "fit_model": "po", "MSE": 0.0808414, "MAD": 0.2305295, "TV": 0.3602959},
    {"n": 100, "sim_data": "ph", "fit_model": "ph", "MSE": 0.08445374, "MAD": 0.2262808, "TV": 0.3614333},
    {"n": 100, "sim_data": "ph", "fit_model": "np", "MSE": 0.08460945, "MAD": 0.2359137, "TV": 0.3431483},
    # n = 100, sim_data = "np"
    {"n": 100, "sim_data": "np", "fit_model": "po", "MSE": 0.4382384, "MAD": 0.5334929, "TV": 0.5966239},
    {"n": 100, "sim_data": "np", "fit_model": "ph", "MSE": 0.4384855, "MAD": 0.5289496, "TV": 0.6166853},
    {"n": 100, "sim_data": "np", "fit_model": "np", "MSE": 0.3271086, "MAD": 0.4542267, "TV": 0.5798078},
    # n = 100, sim_data = "real"
    {"n": 100, "sim_data": "real", "fit_model": "po", "MSE": 1.540615, "MAD": 1.055239, "TV": 0.709665},
    {"n": 100, "sim_data": "real", "fit_model": "ph", "MSE": 1.691484, "MAD": 1.149635, "TV": 0.7318824},
    {"n": 100, "sim_data": "real", "fit_model": "np", "MSE": 1.239083, "MAD": 0.9760289, "TV": 0.6679919},
    
    # n = 500, sim_data = "po"
    {"n": 500, "sim_data": "po", "fit_model": "po", "MSE": 0.03795464, "MAD": 0.1604402, "TV": 0.2513428},
    {"n": 500, "sim_data": "po", "fit_model": "ph", "MSE": 0.03553334, "MAD": 0.1553665, "TV": 0.2456608},
    {"n": 500, "sim_data": "po", "fit_model": "np", "MSE": 0.03621744, "MAD": 0.153854, "TV": 0.2542119},
    # n = 500, sim_data = "ph"
    {"n": 500, "sim_data": "ph", "fit_model": "po", "MSE": 0.01627649, "MAD": 0.1026499, "TV": 0.1641359},
    {"n": 500, "sim_data": "ph", "fit_model": "ph", "MSE": 0.0167663, "MAD": 0.1042484, "TV": 0.1667657},
    {"n": 500, "sim_data": "ph", "fit_model": "np", "MSE": 0.01981645, "MAD": 0.1124745, "TV": 0.1694499},
    # n = 500, sim_data = "np"
    {"n": 500, "sim_data": "np", "fit_model": "po", "MSE": 0.05969665, "MAD": 0.2001962, "TV": 0.3578482},
    {"n": 500, "sim_data": "np", "fit_model": "ph", "MSE": 0.07294462, "MAD": 0.21594, "TV": 0.3634039},
    {"n": 500, "sim_data": "np", "fit_model": "np", "MSE": 0.06601938, "MAD": 0.2027882, "TV": 0.3569291},
    # n = 500, sim_data = "real"
    {"n": 500, "sim_data": "real", "fit_model": "po", "MSE": 1.24663, "MAD": 1.079048, "TV": 0.5771201},
    {"n": 500, "sim_data": "real", "fit_model": "ph", "MSE": 1.284086, "MAD": 1.097552, "TV": 0.5897076},
    {"n": 500, "sim_data": "real", "fit_model": "np", "MSE": 1.112792, "MAD": 1.031674, "TV": 0.46374},
    
    # n = 1000, sim_data = "po"
    {"n": 1000, "sim_data": "po", "fit_model": "po", "MSE": 0.01485118, "MAD": 0.09911655, "TV": 0.1925733},
    {"n": 1000, "sim_data": "po", "fit_model": "ph", "MSE": 0.01959501, "MAD": 0.1096165, "TV": 0.1958018},
    {"n": 1000, "sim_data": "po", "fit_model": "np", "MSE": 0.02023656, "MAD": 0.1146229, "TV": 0.1941967},
    # n = 1000, sim_data = "ph"
    {"n": 1000, "sim_data": "ph", "fit_model": "po", "MSE": 0.007778821, "MAD": 0.06984608, "TV": 0.1250207},
    {"n": 1000, "sim_data": "ph", "fit_model": "ph", "MSE": 0.007064, "MAD": 0.06825603, "TV": 0.1211343},
    {"n": 1000, "sim_data": "ph", "fit_model": "np", "MSE": 0.008162473, "MAD": 0.07145045, "TV": 0.1249989},
    # n = 1000, sim_data = "np"
    {"n": 1000, "sim_data": "np", "fit_model": "po", "MSE": 0.03326955, "MAD": 0.1448028, "TV": 0.283965},
    {"n": 1000, "sim_data": "np", "fit_model": "ph", "MSE": 0.04125971, "MAD": 0.1588127, "TV": 0.28271},
    {"n": 1000, "sim_data": "np", "fit_model": "np", "MSE": 0.02950231, "MAD": 0.1393802, "TV": 0.278772},
    # n = 1000, sim_data = "real"
    {"n": 1000, "sim_data": "real", "fit_model": "po", "MSE": 1.128931, "MAD": 1.047526, "TV": 0.4584095},
    {"n": 1000, "sim_data": "real", "fit_model": "ph", "MSE": 1.153706, "MAD": 1.056565, "TV": 0.4731752},
    {"n": 1000, "sim_data": "real", "fit_model": "np", "MSE": 1.120246, "MAD": 1.047892, "TV": 0.4262749},
]

df = pd.DataFrame(data)

# Create separate figures for each error metric: MSE, MAD, and Average TV distance
metrics = ['MSE', 'MAD', 'TV']
titles = {'MSE': 'Mean Squared Error', 'MAD': 'Mean Absolute Deviation', 'TV': 'Average TV Distance'}

# For each metric, create a 2x2 subplot (one per simulation data type)
for metric in metrics:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axs = axs.flatten()
    sim_types = sorted(df['sim_data'].unique())
    for i, sim in enumerate(sim_types):
        ax = axs[i]
        df_sim = df[df['sim_data'] == sim]
        for fit in sorted(df_sim['fit_model'].unique()):
            df_plot = df_sim[df_sim['fit_model'] == fit].sort_values('n')
            ax.plot(df_plot['n'], df_plot[metric], marker='o', label=f'Fit: {fit}')
        ax.set_title(f'Simulated Data: {sim}')
        ax.set_xlabel('Sample Size (n)')
        ax.set_ylabel(metric)
        ax.legend()
    plt.suptitle(titles[metric])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

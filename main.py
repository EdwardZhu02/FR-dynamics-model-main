# main.py

from utils.data_processing import load_data, process_data
from models.photosynthesis import calculate_net_assimilation
from models.allocation import carbon_allocation, nitrogen_allocation
from models.n_uptake import calculate_n_uptake, optimal_n_rooting_depth
from models.water_uptake import calculate_water_uptake, optimal_water_rooting_depth
from models.n_balance import daily_n_balance
from models.water_balance import daily_water_balance
from models.root_distribution import calculate_root_distribution, root_respiration
from utils.visualization import plot_root_distribution, plot_balance_components


def main():
    # Load and preprocess data
    data = load_data('data/input/')
    processed_data = process_data(data)

    # Calculate net assimilation
    net_assimilation = calculate_net_assimilation(processed_data)

    # Carbon and nitrogen allocation
    carbon_content = carbon_allocation(net_assimilation)
    nitrogen_content = nitrogen_allocation(processed_data)

    # Nitrogen uptake and optimal rooting depth
    n_uptake = calculate_n_uptake(nitrogen_content, processed_data)
    n_optimal_depth = optimal_n_rooting_depth(n_uptake)

    # Water uptake and optimal rooting depth
    water_uptake = calculate_water_uptake(processed_data)
    water_optimal_depth = optimal_water_rooting_depth(water_uptake)

    # Daily nitrogen and water balance
    n_balance = daily_n_balance(n_optimal_depth, processed_data)
    water_balance = daily_water_balance(water_optimal_depth, processed_data)

    # Fine root distribution and respiration
    root_distribution = calculate_root_distribution(n_optimal_depth, water_optimal_depth)
    root_respiration = root_respiration(root_distribution)

    # Visualize results
    plot_root_distribution(root_distribution)
    plot_balance_components(n_balance, water_balance)

    # Save outputs
    # ... code to save outputs ...


if __name__ == '__main__':
    main()

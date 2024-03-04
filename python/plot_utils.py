import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="ticks", palette='muted')

from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler


def plot_training_history(train_losses, test_losses, title='Training Loss over Epochs'):
    """
    Plot the training loss over epochs.

    Parameters:
    - train_losses (list): List containing training loss values for each epoch.
    - title (str): Title for the plot (default is 'Training Loss over Epochs').
    """
    plt.plot(train_losses, label='Training Loss')
    plt.plot(test_losses, label='Test Loss')
    plt.title(title)
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_val(test_labels, test_predictions, rmse_limit=2.8, axes=[7, 16]):
    plt.xlabel('Observed values')
    plt.ylabel('Predicted values')
    plt.xlim(axes)
    plt.ylim(axes)

    mask_color = np.logical_or((test_predictions >= (test_labels + rmse_limit)),
                               (test_predictions <= (test_labels - rmse_limit)))
    per_in = 1 - np.sum(mask_color) / mask_color.size
    color = np.where(mask_color, 'out', f'in ({per_in * 100:.2f}%)')
    sns.scatterplot(x=test_labels.flatten(), y=test_predictions.flatten(),
                    hue=color.flatten(), palette=['#212c3d', '#c05e31'])
    mini = min(test_labels.flatten().min(), test_predictions.flatten().min())
    maxi = max(test_labels.flatten().max(), test_predictions.flatten().max())

    plt.plot([mini, maxi], [mini, maxi], linestyle='--', color='gray')
    plt.plot([mini, maxi], [mini - rmse_limit, maxi - rmse_limit], linestyle='--', color='black')
    plt.plot([mini, maxi], [mini + rmse_limit, maxi + rmse_limit], linestyle='--', color='black')

    # Custom legend labels
    legend_labels = {'out': f'Beyond limit ({np.sum(mask_color)})',
                     'in': f'Inside interval ({per_in * 100:.0f}%)'}
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#c05e31',
                                  markersize=10, label=legend_labels['out']),
                       plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#212c3d',
                                  markersize=10, label=legend_labels['in'])]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.show()


def plot_err(test_labels, test_predictions, bins=35, deviation=1, n_count=50):
    plt.figure()
    error = test_predictions - test_labels
    plt.hist(error, bins)
    plt.xlabel("Prediction Error")
    _ = plt.ylabel("Count")
    plt.xlim([-deviation, deviation])
    plt.ylim([0, n_count])
    plt.show()


def plot_history(history, rmsep, r2, usedForPaper=False):
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    if not usedForPaper:
        plt.title(f'Train/val loss - RMSEP: {rmsep:.3f}; R2: {r2:.3f}')
    plt.ylabel('Loss')
    plt.xlabel('Epochs')
    plt.legend(['train', 'validation'], loc='upper right')
    plt.show()


def plot_prediction(test_labels, test_predictions):
    plt.figure()
    plt.scatter(test_labels, test_predictions)
    # plt.axis('equal')
    plt.legend()
    plt.xlim([0, 80])
    plt.ylim([0, 80])

    plt.xlabel('Observed values')
    plt.ylabel('Predicted values')
    _ = plt.plot([-100, 100], [-100, 100], linestyle='--', color='gray')

    plt.figure()
    error = test_predictions - test_labels
    plt.hist(error, bins=50)
    plt.xlabel("Prediction Error")
    _ = plt.ylabel("Count")


def explore_data(x, variance=.99):
    """
    PCA compute to explore given data
    x: np.array -> the input data to explore
    variance: float -> choose between 95% or 99% on a scale of 0. to 1.
    """
    scaler = MinMaxScaler()
    data_rescaled = scaler.fit_transform(x)

    pca = PCA(n_components=variance)
    X_pca = pca.fit_transform(data_rescaled)
    sns.scatterplot(X_pca)


def plot_df(values, xlabel, fontsize=30, nbr_of_ticks=25):
    df = pd.DataFrame(values)
    df.columns = xlabel.T
    _, ax = plt.subplots(1, figsize=(18,10))

    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(fontsize-7)

    ax.plot(df.transpose())
    plt.xlabel('Wavenumbers cm-1', fontsize=fontsize)
    plt.ylabel("Absorbance", fontsize=fontsize)
    plt.show()

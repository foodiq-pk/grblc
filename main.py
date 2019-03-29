from data_tier_placeholder.image import Image
from matplotlib import pyplot as plt


def plot_src_flux(image_list: [Image]):
    time = []
    src_flux = []

    for img in image_list:
        time.append(img.get_time_jd())
        src_flux.append(img.get_src_flux()[0])
    plt.plot(time, src_flux, ".")


def plot_sky(image_list: [Image]):
    time = []
    sky = []
    for img in image_list:
        time.append(img.get_time_jd())
        sky.append(img.get_sky()[0])
    plt.plot(time, sky, ".")

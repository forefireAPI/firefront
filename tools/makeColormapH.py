import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

def generate_colormap_header(colormaps, num_samples=10, scalemap=[0, 0.9, 1], alpha=[0.3, 0.6, 0.8, 1], cbarPath=None, cbarTextRange=[2, 90]):
    header_content = """
#ifndef COLORMAP_H
#define COLORMAP_H

#include <array>
#include <map>
#include <string>

typedef std::array<unsigned char, 4> ColorEntry;
constexpr size_t colormapSize = {};

""".format(num_samples)

    scale_input = np.linspace(0, 1, len(scalemap))
    alpha_input = np.linspace(0, 1, len(alpha))

    for cmap_name in colormaps:
        try:
            cmap = plt.get_cmap(cmap_name)
        except ValueError:
            continue

        header_content += "static const ColorEntry {}Colormap[colormapSize] = {{\n".format(cmap_name)
        for i in range(num_samples):
            idx = i / (num_samples - 1)
            idx_scaled = np.interp(idx, scale_input, scalemap)
            color = cmap(idx_scaled)
            a_interp = np.interp(idx, alpha_input, alpha)
            r, g, b = [int(255 * x) for x in color[:3]]
            a = int(255 * a_interp)
            header_content += "    {{{}, {}, {}, {}}},\n".format(r, g, b, a)
        header_content += "};\n\n"

        if cbarPath is not None:
            fig, ax = plt.subplots(figsize=(2, 10))
            bounds = np.linspace(0, 1, num_samples + 1)
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
            cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, ticks=bounds)
            scaled_positions = [np.interp(i / num_samples, scale_input, scalemap) * (cbarTextRange[1] - cbarTextRange[0]) + cbarTextRange[0] for i in range(num_samples + 1)]
            cb.set_ticklabels([f"{l:.2f}" for l in scaled_positions])
            plt.tight_layout()
            plt.savefig(f"{cbarPath}{cmap_name}.png", bbox_inches='tight', dpi=300)
            plt.close()

    header_content += "static const std::map<std::string, const ColorEntry*> colormapMap = {\n"
    for cmap_name in colormaps:
        header_content += '    {{"{}", {}Colormap}},\n'.format(cmap_name, cmap_name)
    header_content += "};\n"
    header_content += "#endif // COLORMAP_H\n"
    return header_content

# Usage
colormaps = ['hot', 'viridis', 'plasma', 'plasma_r', 'coolwarm', 'grey', 'hot_r', 'spring', 'jet', 'turbo']
header_content = generate_colormap_header(colormaps, 50, scalemap=[0, 0.2, 0.5, 1], alpha=[0.3, 0.5, 0.9, 1], cbarPath="cmaps/", cbarTextRange=[0, 8000])
with open("../src/colormap.h", "w") as file:
    file.write(header_content)

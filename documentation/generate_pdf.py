#!/usr/bin/env python3
"""Generate figures and compile the SAFB documentation PDF.

Generates:
1. architecture_diagram.png - Python to C module mapping
2. star_generation.png - Star generation visualization
3. field_morphology_gyroid.png - Gyroid field from Ia-3d
4. field_morphology_bcc.png - BCC-like field from Pm-3m

Then compiles the LaTeX document.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from fractions import Fraction

DOC_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = DOC_DIR

def generate_architecture_diagram():
    """Generate Python-to-C module mapping diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    ax.set_title('SAFB Library Architecture: Python → C Mapping', fontsize=14, fontweight='bold', pad=20)

    # Python modules (left)
    py_modules = [
        ('domain.py', 'Data structures'),
        ('symmetry_ops.py', 'Symmetry ops'),
        ('space_group_plane_family.py', 'Basis construction'),
        ('initializers.py', 'Amplitude assignment'),
        ('field.py', 'Field generation'),
        ('engine.py', 'High-level API'),
        ('Analytic.py', 'Analytical calc.'),
    ]

    # C modules (right)
    c_modules = [
        ('include/domain.h\nsrc/domain.c', 'Data structures'),
        ('include/symmetry_ops.h\nsrc/symmetry_ops.c', 'Symmetry ops'),
        ('include/basis.h\nsrc/basis.c', 'Basis construction'),
        ('include/initializers.h\nsrc/initializers.c', 'Amplitude assignment'),
        ('include/field.h\nsrc/field.c', 'Field generation'),
        ('include/engine.h\nsrc/engine.c', 'High-level API'),
        ('include/analytic.h\nsrc/analytic.c', 'Analytical calc.'),
    ]

    # Draw Python modules
    for i, (name, desc) in enumerate(py_modules):
        y = 8.5 - i * 1.2
        box = FancyBboxPatch((0.5, y - 0.35), 3.5, 0.7,
                             boxstyle="round,pad=0.05",
                             facecolor='#E3F2FD', edgecolor='#1565C0', linewidth=2)
        ax.add_patch(box)
        ax.text(2.25, y, name, ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(2.25, y - 0.25, desc, ha='center', va='center', fontsize=7, style='italic')

    # Draw C modules
    for i, (name, desc) in enumerate(c_modules):
        y = 8.5 - i * 1.2
        box = FancyBboxPatch((6, y - 0.35), 3.5, 0.7,
                             boxstyle="round,pad=0.05",
                             facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2)
        ax.add_patch(box)
        ax.text(7.75, y, name, ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(7.75, y - 0.25, desc, ha='center', va='center', fontsize=7, style='italic')

    # Draw arrows
    for i in range(len(py_modules)):
        y = 8.5 - i * 1.2
        arrow = FancyArrowPatch((4.2, y), (5.8, y),
                               arrowstyle='->', mutation_scale=20,
                               color='#666666', linewidth=1.5)
        ax.add_patch(arrow)

    # Title labels
    ax.text(2.25, 9.5, 'Python Source\n(/sandbox/Sg_init)', ha='center', fontsize=10, fontweight='bold')
    ax.text(7.75, 9.5, 'C Target\n(/sandbox/hermes_SAFB_migration)', ha='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'architecture_diagram.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Generated: architecture_diagram.png")


def generate_star_generation():
    """Generate star generation visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Example: 4-fold star from (1,0,0) in cubic lattice
    # Point group rotations (4-fold)
    rots_2d = [
        [[1, 0], [0, 1]],
        [[0, -1], [1, 0]],
        [[-1, 0], [0, -1]],
        [[0, 1], [-1, 0]],
    ]
    seed = np.array([1, 0])

    for ax_idx, (ax, seed_label) in enumerate(zip(axes.flat[:4], [(1,0), (1,1), (2,0), (1,0,-1)])):
        ax.set_aspect('equal')
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)

        if len(seed_label) == 2:
            seed_arr = np.array(seed_label)
        else:
            seed_arr = np.array(seed_label[:2])

        # Apply rotations
        vectors = []
        for R in rots_2d:
            R_arr = np.array(R)
            v = R_arr @ seed_arr
            vectors.append(v)
        # Add negatives
        for v in vectors:
            vectors.append(-v)

        vectors = np.array(vectors)

        # Draw vectors from origin
        for i, v in enumerate(vectors):
            ax.arrow(0, 0, v[0], v[1], head_width=0.15, head_length=0.1,
                    fc='steelblue', ec='steelblue', alpha=0.7, linewidth=1.5)

        # Draw circle at radius |q|
        radius = np.linalg.norm(seed_arr)
        circle = plt.Circle((0, 0), radius, fill=False, color='red', linestyle='--', linewidth=1.5)
        ax.add_patch(circle)

        # Label seed
        ax.plot(seed_arr[0], seed_arr[1], 'ro', markersize=10, label='Seed')
        ax.plot(0, 0, 'ko', markersize=6)

        ax.set_title(f'Star({seed_label})\nMultiplicity: {len(vectors)}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)

        if ax_idx == 0:
            ax.legend(loc='upper right', fontsize=8)

    # Remove unused subplot
    axes[1, 1].set_visible(False)

    fig.suptitle('Star Generation: Symmetry-Equivalent Vectors', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'star_generation.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Generated: star_generation.png")


def generate_field_morphology_gyroid():
    """Generate gyroid field from Ia-3d (simulated)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Simulate gyroid-like field using Fourier modes
    # Gyroid: ψ = 2[cos(q₁·r)cos(q₂·r) + cos(q₂·r)cos(q₃·r) + cos(q₃·r)cos(q₁·r)]
    # where q₁, q₂, q₃ are the three cubic face diagonals

    Nx, Ny, Nz = 64, 64, 64
    x = np.linspace(0, 2*np.pi, Nx)
    y = np.linspace(0, 2*np.pi, Ny)
    z = np.linspace(0, 2*np.pi, Nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Gyroid equation (simplified)
    q1 = np.array([1, 1, 0])
    q2 = np.array([0, 1, 1])
    q3 = np.array([1, 0, 1])

    psi = (np.cos(q1[0]*X + q1[1]*Y) * np.cos(q2[0]*Y + q2[1]*Z) +
           np.cos(q2[0]*Y + q2[1]*Z) * np.cos(q3[0]*X + q3[1]*Z) +
           np.cos(q3[0]*X + q3[1]*Z) * np.cos(q1[0]*X + q1[1]*Y))

    # Normalize to [0, 1]
    psi = (psi - psi.min()) / (psi.max() - psi.min())

    # Plot slices
    mid = Nz // 2
    axes[0, 0].imshow(psi[:, :, mid], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[0, 0].set_title('Gyroid ψ(x,y,z=π)', fontsize=10, fontweight='bold')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')

    axes[0, 1].imshow(psi[:, mid, :], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[0, 1].set_title('Gyroid ψ(x,z=π,y)', fontsize=10, fontweight='bold')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Z')

    axes[1, 0].imshow(psi[mid, :, :], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[1, 0].set_title('Gyroid ψ(y=π,x,z)', fontsize=10, fontweight='bold')
    axes[1, 0].set_xlabel('Y')
    axes[1, 0].set_ylabel('Z')

    # Histogram
    axes[1, 1].hist(psi.ravel(), bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    axes[1, 1].set_title('Field Distribution', fontsize=10, fontweight='bold')
    axes[1, 1].set_xlabel('ψ value')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].axvline(psi.mean(), color='red', linestyle='--', label=f'Mean={psi.mean():.3f}')
    axes[1, 1].legend()

    fig.suptitle('SAFB Gyroid Field (Ia-3d) — Simulated Morphology', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'field_morphology_gyroid.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Generated: field_morphology_gyroid.png")


def generate_field_morphology_bcc():
    """Generate BCC-like field from Pm-3m (simulated)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    Nx, Ny, Nz = 64, 64, 64
    x = np.linspace(0, 2*np.pi, Nx)
    y = np.linspace(0, 2*np.pi, Ny)
    z = np.linspace(0, 2*np.pi, Nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # BCC-like field: sum of cosines with cubic symmetry
    psi = (np.cos(X) + np.cos(Y) + np.cos(Z) +
           0.5 * (np.cos(2*X) + np.cos(2*Y) + np.cos(2*Z)))

    # Normalize to [0, 1]
    psi = (psi - psi.min()) / (psi.max() - psi.min())

    # Plot slices
    mid = Nz // 2
    axes[0, 0].imshow(psi[:, :, mid], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[0, 0].set_title('BCC ψ(x,y,z=π)', fontsize=10, fontweight='bold')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')

    axes[0, 1].imshow(psi[:, mid, :], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[0, 1].set_title('BCC ψ(x,z=π,y)', fontsize=10, fontweight='bold')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Z')

    axes[1, 0].imshow(psi[mid, :, :], extent=[0, 2*np.pi, 0, 2*np.pi],
                     cmap='RdBu_r', origin='lower', aspect='equal')
    axes[1, 0].set_title('BCC ψ(y=π,x,z)', fontsize=10, fontweight='bold')
    axes[1, 0].set_xlabel('Y')
    axes[1, 0].set_ylabel('Z')

    # Histogram
    axes[1, 1].hist(psi.ravel(), bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    axes[1, 1].set_title('Field Distribution', fontsize=10, fontweight='bold')
    axes[1, 1].set_xlabel('ψ value')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].axvline(psi.mean(), color='red', linestyle='--', label=f'Mean={psi.mean():.3f}')
    axes[1, 1].legend()

    fig.suptitle('SAFB BCC Field (Pm-3m) — Simulated Morphology', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'field_morphology_bcc.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Generated: field_morphology_bcc.png")


def compile_pdf():
    """Compile the LaTeX document to PDF."""
    import subprocess

    os.chdir(DOC_DIR)

    # Install texlive if needed
    result = subprocess.run(['which', 'pdflatex'], capture_output=True)
    if result.returncode != 0:
        print("Installing texlive...")
        subprocess.run(['conda', 'install', '-c', 'conda-forge',
                       'texlive-core', 'texlive-latexextra', 'texlive-science', '-y'],
                      env={**os.environ, 'CONDA_NO_PLUGINS': 'true',
                           'CONDA_OVERRIDE_CUDA': '', 'solver': 'classic'},
                      check=True)

    # First pass
    result = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'main.tex'],
                           capture_output=True, text=True)
    if result.returncode != 0:
        print("First pass warnings/errors:")
        print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)

    # Second pass (for references)
    result = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'main.tex'],
                           capture_output=True, text=True)
    if result.returncode != 0:
        print("Second pass warnings/errors:")
        print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)

    # Check if PDF was created
    pdf_path = os.path.join(DOC_DIR, 'main.pdf')
    if os.path.exists(pdf_path):
        print(f"PDF generated: {pdf_path}")
        return True
    else:
        print("PDF generation failed. Check LaTeX output above.")
        return False


def main():
    print("=" * 60)
    print("SAFB Documentation: Figure Generation & PDF Compilation")
    print("=" * 60)

    # Generate figures
    print("\n--- Generating Figures ---")
    generate_architecture_diagram()
    generate_star_generation()
    generate_field_morphology_gyroid()
    generate_field_morphology_bcc()

    # Compile PDF
    print("\n--- Compiling PDF ---")
    compile_pdf()

    print("\nDone!")


if __name__ == "__main__":
    main()

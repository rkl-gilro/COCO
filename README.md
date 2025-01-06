# Colour Constancy in Virtual Reality

The two indoor and outdoor setups in the study focus on investigating colour constancy mechanisms in virtual reality environments. 

### Indoor Setup
- **Environment:** A realistic office-like space with objects such as a desk, chair, bookshelf, plants, and decorative elements.
- **Illumination:** Single point light source positioned above the observer’s VR actor in the centre of the room.
- **Scene Interactions:** Observers adapted to five illuminants (neutral, blue, green, yellow, red).
- **Key Experiment Aspects:**
    1. The achromatic reference object (a lizard) was placed on a painting under neutral illumination.
    2. Competing objects of varying colours were introduced under different illuminants.


### Outdoor Setup
- **Environment:** A naturalistic forest scene with features like trees, cliffs, a lake, rocks, moss, grass, and flowers.
- **Illumination:** A directional light paired with a skylight to simulate natural conditions, with adjustments based on illuminant colours.
- **Scene Interactions:** Observers adapted to five illuminants (neutral, blue, green, yellow, red).
- **Key Experiment Aspects:**
    1. A reference lizard was presented on moss near the cliff under neutral illumination.
    2. Competing lizards were placed around the environment in random locations.

Both setups utilized VR to maintain experimental control, allowing for precise manipulations of illumination and visual cues while preserving immersive, realistic conditions. The study aimed to understand how mechanisms like local surround, maximum flux, and spatial mean influence color constancy under varied lighting.

## Scene with one light source (indoor and outdoor)

### Local Surround mechanism


<p align="center" width="100%">
    <img width="63%" src="/pictures/localsurround1.png">
    <img width="33%" src="/pictures/localsurround2.png">
</p>

The local surround mechanism examines the influence of immediate surroundings on colour constancy.

**Purpose:**

To test the impact of silencing the local surround cue, which typically stabilizes colour perception through contrasts between the object and its surrounding background.

**Implementation:**

The target object (a lizard) was placed on a consistent, self-illuminated leaf with a fixed chromaticity regardless of illuminant colour.

**Effect:**

This manipulation ensured that the immediate background contrast (local surround) remained constant and neutral, isolating the role of this cue in colour constancy.



## Scene with two light sources (indoor illumination + sunlight coming through window)

<p align="center" width="100%">
    <img width="63%" src="/pictures/sceneCC2.png">
</p>

**Objective**
The study investigates colour constancy under conditions with two illuminants, focusing on tasks like selection, asymmetric matching, and achromatic adjustment.

**Methodology**
  1. Selection Task:

      ` Participants identify a "Lonely Kula" object under two neutral illuminants and later in a coloured room using a controller.
      Objects include cubes and spheres.`

  2. Asymmetric Matching:

        ` Adjust the test object’s colour to match a reference object by sliding a controller’s trackpad.`

  3. Achromatic Adjustment:

        `Adjust the test object’s colour to appear achromatic under both neutral and coloured illuminants.
        Two object positions were tested.`

**Quantification**

Colour Constancy Index (CCI): Measures performance by comparing the differences between test and reference matches under neutral and coloured illuminants.

**Preliminary Results**

  1. Selection Task:
      - Achieved the highest CCI scores.
      - Shows preference for cubes over spheres.

  2. Asymmetric Matching:
     - Challenging for participants.
     - Poor CCI due to difficulty and coarse adjustments.
 
  3. Achromatic Adjustment:
     - No significant differences between positions.
     - Moderate CCI performance.

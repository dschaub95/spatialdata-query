// this script is broken as it leads to a out-of-memory error
import { ParquetLoader } from '@loaders.gl/parquet';
import { load } from '@loaders.gl/core';
import { promises as fs } from 'fs';

async function processParquetFile() {
    try {
        // Read the Parquet file into a buffer
        const fileBuffer = await fs.readFile('random_points.parquet');

        // Load the Parquet data from the buffer
        const data = await load(fileBuffer, ParquetLoader);

        // Extract (x, y) coordinates into an array
        const coordinates = data.map(row => [row.x, row.y]);

        // Print the shape of the coordinates array
        console.log(`Shape of coordinates array: [${coordinates.length}, 2]`);
    } catch (error) {
        console.error('Error loading Parquet file:', error);
    }
}

processParquetFile();

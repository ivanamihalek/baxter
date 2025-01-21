import { unstable_cache } from 'next/cache';

// Define the type for the single JSON entry
interface FingerprintAndDecoys {
    fingerprint: string;
    decoy1: string;
    decoy2: string;
}

// Define the type for the data returned by fetchData
type DataType = FingerprintAndDecoys; // Change to a single object

// Define an interface for the cache configuration
interface CacheConfig<DataType> {
    fetchData: (url: string) => Promise<DataType>; // Function type for fetching data
    keyParts: string[];                             // Array of strings for key parts
    options: CacheOptions;                          // Options using the CacheOptions interface
}

// Define an interface for the options parameter
interface CacheOptions {
    revalidate?: number; // Optional revalidation interval in seconds
    tags?: string[];     // Optional array of tags for cache invalidation
}

async function fetchData(url: string): Promise<DataType> {
 
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error('Network response was not ok.');
    }

    return response.json();
}

// Create a cached version of fetchData using unstable_cache with the defined interface
const cacheConfig: CacheConfig<DataType> = {
    fetchData,
    keyParts: ['external-data'],
    options: {
        revalidate: 2, // Cache for two seconds or any time you like
        tags: ['external-data'],  // Tag for invalidation
    },
};

const cachedFetchData = unstable_cache(
    cacheConfig.fetchData,
    cacheConfig.keyParts,
    cacheConfig.options
);

export default async function Home() {
    let fpd: FingerprintAndDecoys | null = null; // Change to a single object or null

    try {
        fpd = await cachedFetchData('http://localhost:8000/bad_bac_exercise/fingerprint/');
        // menuItem = await cachedFetchData('http://baxter-backend-ctnr:8000/bad_bac_exercise/arm/');
    } catch (err) {
        console.error(err);
    }

    // Handle case where menuItem is null
    if (fpd === null) {
        return <div>No data available.</div>; // Or handle it as you see fit
    }

    const seqArray = [fpd.fingerprint, fpd.decoy1, fpd.decoy2]
    const shuffledSeqArray = seqArray.sort(() => Math.random() - 0.5);
    return (
        <div>
            <h1>Sequence snippets</h1>
            <ul>
                    {shuffledSeqArray.map((item, index) => (
                    <li key={index}>{item}</li> // Use index as key for simplicity here
                    ))}
            </ul>
        </div>
    );
}

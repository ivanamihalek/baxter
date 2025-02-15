/**
 * @typedef {object} FingerprintAndDecoys
 * @property {string} fingerprint
 * @property {string} decoy1
 * @property {string} decoy2
 */

/**
 * @typedef {object} CacheConfig
 * @property {function(string): Promise<FingerprintAndDecoys>} fetchData - Function to fetch data
 * @property {string[]} keyParts - Array of strings for key parts
 * @property {CacheOptions} options - Options for cache configuration
 */

/**
 * @typedef {object} CacheOptions
 * @property {number} [revalidate] - Optional revalidation interval in seconds
 * @property {string[]} [tags] - Optional array of tags for cache invalidation
 */

import {unstable_cache} from "next/cache";

/**
 * @async
 * @function fetchData
 * @param {string} url
 * @returns {Promise<FingerprintAndDecoys>}
 */
async function fetchData(url) {
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error('Network response was not ok.');
    }
    return response.json();
}

/**
 * @type {CacheConfig}
 */
const cacheConfig = {
    fetchData: fetchData,
    keyParts: ['external-data'],
    options: {
        revalidate: 2, // Cache for two seconds
        tags: ['external-data'],  // Tag for invalidation
    },
};

const cachedFetchData = unstable_cache(
    cacheConfig.fetchData,
    cacheConfig.keyParts,
    cacheConfig.options
);

/**
 * @async
 * @function Home
 * @returns {JSX.Element}
 */
export default async function Home() {
    /** @type {FingerprintAndDecoys | null} */
    let fpd = null;

    try {
        // todo - set this address somehow so it's different in production and dev
        fpd = await cachedFetchData('http://baxter-backend-ctnr:8000/bad_bac_exercise/fingerprint');
    } catch (err) {
        console.error(err);
    }

    if (fpd === null) {
        return <div>No data available.</div>;
    }

    const seqArray = [fpd.fingerprint, fpd.decoy1, fpd.decoy2];
    const shuffledSeqArray = seqArray.sort(() => Math.random() - 0.5);

    return (
        <div>
            <h1>Sequence snippets</h1>
            <ul>
                {shuffledSeqArray.map((item, index) => (
                    <li key={index}>{item}</li>
                ))}
            </ul>
        </div>
    );
}

export var backend_url = "";
if (process.env.NODE_ENV === 'production') {
    var backend_url = process.env.REACT_APP_API_URL;
}


export function getHDNeighbors(selectedPoints, kNeighbors, metric) {
    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 0) {
            fetch(`${backend_url}/backend/computeHDNeighbors`, {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/computeHDNeighbors got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                        resolve(response);
                    } else {
                        response.json().then(data => {
                            let newOpacities = data["binary"]
                            resolve(newOpacities);
                        })
                    }
                })
        }
        else {
            resolve([]);
        }
    })
}


export function explainCluster(selectedPoints, explainabilityMethod) {
    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 0) {
            fetch(`${backend_url}/backend/explainCluster`, {
                method: "POST",
                body: JSON.stringify({
                    points: selectedPoints,
                    selection_name: explainabilityMethod,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/explainCluster got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                        resolve({"features": [], "higher_mean": []});
                    } else {
                        response.json().then(data => {
                            resolve(data);
                        })
                    }
                })
        }
        else {
            resolve({"features": [], "higher_mean": []});
        }
    })
}

export function compareClusters(selectionA, selectionB) {
    return new Promise((resolve, reject) => {
        if (selectionA.length > 0 && selectionB.length > 0) {
            fetch(`${backend_url}/backend/compareClusters`, {
                method: "POST",
                body: JSON.stringify({
                    selectionA: selectionA,
                    selectionB: selectionB,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/compareClusters got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                        resolve([]);
                    } else {
                        response.json().then(data => {
                            resolve(data);
                        })
                    }
                })
        }
        else {
            resolve([]);
        }
    })
}

export function getFeatureValues(featureName) {
    return new Promise((resolve, reject) => {
        fetch(`${backend_url}/backend/featureValues/${featureName}`,
            {
                method: "GET",
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
            .then(response => {
                if (!response.ok) {
                    console.log(`/backend/featureValues/${featureName} got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                    resolve([]);
                } else {
                    response.json().then(data => {
                        resolve(data["result"]);
                    })
                }
            })
        })
    }